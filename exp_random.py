#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Дополнительная серия экспериментов (п. 3.1 статьи): 436 пар случайных автоматов.

Для каждой пары (A, B) с N = K состояниями и |I| = |U| = |O| = p строятся
  MB1 - обход графа переходов Spec = A∙B + различающие последовательности для мутантов
        первого порядка B с ошибками переходов (новая последовательность после сброса);
  MBC - покрытие переходов B, используемых в композиции, + различающие последовательности
        для тех же мутантов (с продолжением существующих последовательностей).
Проверяется полнота обоих тестов относительно всех обнаружимых мутантов первого порядка
(ошибки переходов и выходов) и выборки до 1000 обнаружимых мутантов второго порядка
(ошибки переходов).

Пары: p = 2, 3, 4 и N = K = 3, 5, ..., 15 - по 20 пар; при p = 2 еще N = K = 20, 30 (по 5 пар)
и N = K = 40, 50 (по 3 пары). Пара k строится от зерна 1000003*N + 7919*p + k
(как в fsm_generator.py).

    python exp_random.py [файл результатов] [-j число процессов]
Время счета около минуты на двух ядрах.
"""
import argparse
import json
import os
import random
import time
from multiprocessing import Pool

from fsmnet import Net, rand_fsm, mb, length, second_order_sample


def grid():
    g = []
    for p in (2, 3, 4):
        for n in range(3, 16, 2):
            g.append((n, n, p, 20))
    g += [(20, 20, 2, 5), (30, 30, 2, 5), (40, 40, 2, 3), (50, 50, 2, 3)]
    return g


def experiment(args):
    N, K, p, seed, so_size = args
    t_start = time.time()
    rng = random.Random(seed)
    A = rand_fsm(N, p, p, rng)
    B = rand_fsm(K, p, p, rng)
    net = Net(A, B, p)
    rec = dict(N=N, K=K, p=p, seed=seed, r=net.r, reachB=len(net.reachB), used=len(net.used))
    RTf = net.first_order(('T',))
    RT = [net.mutant(f) for f in RTf]
    # классификация необнаружимых мутантов первого порядка с ошибками переходов
    unreachB = unusedU = masked = 0
    for f, M in zip(RTf, RT):
        b, u = f[0][0], f[0][1]
        if (b, u) not in net.used:
            if b not in net.reachB: unreachB += 1
            else: unusedU += 1
        elif not net.distinguishable(M):
            masked += 1
    rec.update(nRT=len(RT), unreachB=unreachB, unusedU=unusedU, masked=masked)
    ALLdet = [M for M in (net.mutant(f) for f in net.first_order(('T', 'O'))) if net.distinguishable(M)]
    rec['n1'] = len(ALLdet)

    tour = net.cpp_tour()
    rec['tour'] = length(tour)
    T1, a1 = mb(net, tour, RT, ext=False)
    rec['MB1'] = length(T1); rec['MB1_added'] = a1
    cov = net.comp_cover()
    rec['cover'] = length(cov)
    TC, ac = mb(net, cov, RT, ext=True)
    rec['MBC'] = length(TC); rec['MBC_added'] = ac
    tests = {'MB1': T1, 'MBC': TC}
    for name, T in tests.items():
        rec[name + '_M1'] = sum(net.detects(T, M) for M in ALLdet) / len(ALLdet) if ALLdet else 1.0
    so = second_order_sample(net, rng, so_size)
    rec['n2'] = len(so)
    for name, T in tests.items():
        rec[name + '_2nd'] = sum(net.detects(T, M) for f, M in so) / len(so) if so else 1.0
    rec['time'] = round(time.time() - t_start, 3)
    return rec


def main():
    ap = argparse.ArgumentParser(description='Дополнительная серия экспериментов: 436 пар случайных автоматов')
    ap.add_argument('out', nargs='?', default=os.path.join('results', 'random_436.json'), help='файл результатов (по умолчанию %(default)s)')
    ap.add_argument('-j', '--procs', type=int, default=os.cpu_count() or 1, help='число процессов (по умолчанию - число ядер)')
    a = ap.parse_args()
    out, procs = a.out, a.procs
    jobs = [(N, K, p, 1000003 * N + 7919 * p + k, 1000) for (N, K, p, npairs) in grid() for k in range(npairs)]
    jobs.sort(key=lambda j: -j[0])          # большие пары первыми - равномернее загрузка процессов
    t0 = time.time()
    with Pool(procs) as pool:
        recs = pool.map(experiment, jobs, chunksize=1)
    recs.sort(key=lambda r: (r['p'], r['N'], r['seed']))
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(recs, f, indent=0)
    print(f'пар: {len(recs)}, время {time.time() - t0:.0f} с, результаты: {out}')


if __name__ == '__main__':
    main()
