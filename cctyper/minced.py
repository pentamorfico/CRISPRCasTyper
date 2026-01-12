import os
import logging
import math
import random

import statistics as st
from diced import scan

from joblib import Parallel, delayed

# Define the CRISPR class
class CRISPR(object):
    count = 0
    def __init__(self, sequence, exact_stats):
        self.sequence = sequence.rstrip()
        CRISPR.count += 1
        self.crispr = '{}_{}'.format(self.sequence, CRISPR.count)
        self.repeats = []
        self.spacers = []
        self.exact = exact_stats
    def setPos(self, start, end):
        self.start = int(start)
        self.end = int(end)
    def addRepeat(self, repeat):
        self.repeats.append(repeat.rstrip())
    def addSpacer(self, spacer):
        put_spacer = spacer.rstrip()
        if len(put_spacer) > 0:
            self.spacers.append(put_spacer)
    def getConsensus(self):
        self.cons = max(set(self.repeats), key = self.repeats.count) 
    def identity(self, i, j, sqlst):
        # Equivalent to pairwise2.globalxx: identity based on LCS length over max length
        a = sqlst[i]
        b = sqlst[j]
        if not a or not b:
            return 0
        m, n = len(a), len(b)
        prev = [0] * (n + 1)
        for ca in a:
            curr = [0]
            for j, cb in enumerate(b, start=1):
                if ca == cb:
                    curr.append(prev[j - 1] + 1)
                else:
                    curr.append(max(prev[j], curr[-1]))
            prev = curr
        lcs_len = prev[-1]
        return (lcs_len / float(max(m, n))) * 100
    def identLoop(self, seqs, threads):
        if self.exact:
            sqr = range(len(seqs))
        else:
            if len(seqs) > 10:
                n_samp = 10
            else:
                n_samp = len(seqs)
            sqr = random.sample(range(len(seqs)), n_samp)
        idents = Parallel(n_jobs=threads)(delayed(self.identity)(k, l, seqs) for k in sqr for l in sqr if k > l)
        return(st.mean(idents))
    def stats(self, threads, rep_id, spa_id, spa_sem):
        if len(self.spacers) > 1:
            self.spacer_identity = round(self.identLoop(self.spacers, threads), 1)
            self.spacer_len = round(st.mean([len(x) for x in self.spacers]), 1)
            self.spacer_sem = round(st.stdev([len(x) for x in self.spacers])/math.sqrt(len(self.spacers)), 1)
        else:
            self.spacer_identity = 0
            self.spacer_len = len(self.spacers[0])
            self.spacer_sem = 0
        self.repeat_identity = round(self.identLoop(self.repeats, threads), 1)
        self.repeat_len = round(st.mean([len(x) for x in self.repeats]), 1)
        self.trusted = (self.repeat_identity > rep_id) & (self.spacer_identity < spa_id) & (self.spacer_sem < spa_sem)

class Minced(object):
    
    def __init__(self, obj):
        self.master = obj
        for key, val in vars(obj).items():
            setattr(self, key, val)

    def run_minced(self) -> None:

        if not self.redo:
            logging.info('Predicting CRISPR arrays with diced')

            random.seed(self.seed)
            crisprs = []

            for contig, seq in self.seq_dict.items():
                sequence_str = str(seq)
                try:
                    crispr_calls = list(scan(sequence_str))
                except Exception as exc:
                    logging.error('diced failed on contig %s: %s', contig, exc)
                    continue

                for crispr_call in crispr_calls:
                    try:
                        crisp_tmp = CRISPR(contig, self.exact_stats)
                        seq_len = len(sequence_str)

                        start = max(1, min(crispr_call.start + 1, seq_len))  # convert to 1-based, clamp
                        end = max(start, min(crispr_call.end, seq_len))
                        crisp_tmp.setPos(start, end)

                        # Materialize repeats/spacers once; diced can panic on lazy access if inconsistent
                        bad_repeats = 0
                        bad_spacers = 0

                        # diced can panic on some indices; step through safely
                        for ridx in range(len(crispr_call.repeats)):
                            try:
                                rep = crispr_call.repeats[ridx]
                            except BaseException:
                                bad_repeats += 1
                                continue
                            rs = max(0, min(rep.start, seq_len))
                            re = max(rs, min(rep.end, seq_len))
                            if re > rs:
                                crisp_tmp.addRepeat(sequence_str[rs:re])

                        for sidx in range(len(crispr_call.spacers)):
                            try:
                                spa = crispr_call.spacers[sidx]
                            except BaseException:
                                bad_spacers += 1
                                continue
                            ss = max(0, min(spa.start, seq_len))
                            se = max(ss, min(spa.end, seq_len))
                            if se > ss:
                                crisp_tmp.addSpacer(sequence_str[ss:se])

                        if bad_repeats or bad_spacers:
                            logging.error(
                                'diced returned %s bad repeats and %s bad spacers on %s; skipped those elements',
                                bad_repeats, bad_spacers, contig
                            )

                        if not crisp_tmp.repeats:
                            logging.warning('diced produced a call with no repeats on %s; skipping', contig)
                            continue

                        crisp_tmp.getConsensus()
                        crisp_tmp.stats(self.threads, self.repeat_id, self.spacer_id, self.spacer_sem)
                        crisprs.append(crisp_tmp)
                    except BaseException as exc:
                        logging.error('failed to process diced call on contig %s: %s', contig, exc)
                        continue

            self.crisprs = crisprs
            
            # Write results
            self.write_crisprs()
            self.write_spacers()

    def write_crisprs(self) -> None:
       
        header = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Contig',
                                                                    'CRISPR',
                                                                    'Start',
                                                                    'End',
                                                                    'Consensus_repeat',
                                                                    'N_repeats',
                                                                    'Repeat_len',
                                                                    'Spacer_len_avg',
                                                                    'Repeat_identity',
                                                                    'Spacer_identity',
                                                                    'Spacer_len_sem',
                                                                    'Trusted')
       
        def write_crisp(handle, cris):
            handle.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(cris.sequence,
                                                   cris.crispr,
                                                   cris.start,
                                                   cris.end,
                                                   cris.cons,
                                                   len(cris.repeats),
                                                   cris.repeat_len,
                                                   cris.spacer_len,
                                                   cris.repeat_identity,
                                                   cris.spacer_identity,
                                                   cris.spacer_sem,
                                                   cris.trusted))
            
        f = open(self.out+'crisprs_all.tab', 'w')
        f.write(header)
        for crisp in self.crisprs:
            write_crisp(f, crisp)
        f.close()

    def write_spacers(self) -> None:
        
        if len(self.crisprs) > 0:
            os.mkdir(self.out+'spacers')
            for crisp in self.crisprs:
                f = open(self.out+'spacers/{}.fa'.format(crisp.crispr), 'w')
                n = 0
                for sq in crisp.spacers:
                    n += 1
                    f.write('>{}:{}\n'.format(crisp.crispr, n))
                    f.write('{}\n'.format(sq))

                f.close()
