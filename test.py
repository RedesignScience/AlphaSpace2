import numpy as np
import multiprocessing as mp
import timeit
import time

def multiprocess(function, argslist, ncpu):
    total = len(argslist)
    done = 0
    result_queue = mp.Queue(0)
    jobs = []
    res = []
    while argslist != []:
        if len(mp.active_children()) < ncpu:
            p = mp.Process(target=function, args=(result_queue, argslist.pop(),))
            jobs.append(p)
            p.start()
            done += 1
            print("\r",float(done)/total*100,"%") #here is to keep track
        # here comes my emptying step
        if len(jobs) == 500:
            tmp = [result_queue.get() for p in jobs]
            for r in tmp:
                res.append(r)
            result_queue = mp.Queue(0)
            jobs = []

    tmp = [result_queue.get() for p in jobs]
    for r in tmp:
        res.append(r)
    return res

def test_func( x):
    time.sleep(3)
    return x**2


def test_wrapper(q,x):
    q.put(test_func(x))



if __name__ == '__main__':
    start = timeit.default_timer()
    x_list = [i for i in range(5)]
    result = multiprocess(test_wrapper,x_list,2)
    print(result)
    mp_end_time = timeit.default_timer()
    print(
        mp_end_time - start
    )
    results = []
    for i in range(5):
        result.append(test_func(i))
    print(timeit.default_timer()- mp_end_time)


