__author__ = 'clyde'


from itertools import islice, chain
from multiprocessing.pool import Pool


def batch(iterable, size):
    sourceiter = iter(iterable)
    while True:
        batchiter = islice(sourceiter, size)
        yield chain([batchiter.next()], batchiter)

def parallel_score_it(chunk_it, score_f):
    p = Pool(2)
    for chunk in chunk_it:
        for weighted_accuracy in p.map(score_f, chunk):
            yield weighted_accuracy
    p.close()