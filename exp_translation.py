#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Трансляция теста, построенного для хвостовой компоненты как независимого автомата (п. 3.2 статьи):
серия из 350 пар для оценки полноты транслированного теста.

Для каждой пары (A, B), N = K, |I| = |U| = |O| = p:
  * тест для B в изоляции: W-тест (граница числа состояний K) при N = K <= 4,
    иначе MB12 автомата B (обход графа переходов B + различающие последовательности
    для мутантов B первого и второго порядка);
  * каждая последовательность транслируется через A (fsmnet.translate): до B доходит
    наибольший префикс, который A может породить из начального состояния;
  * доля последовательностей, транслируемых целиком, и доля входных символов, доходящих до B;
  * полнота транслированного теста и теста MB1 (построенного по Spec = A∙B) относительно
    обнаружимых мутантов B первого порядка (все) и второго и третьего порядков (выборки до 1000)
    с ошибками переходов.

Пары: p = 2, N = K = 2..15; p = 3, N = K = 3..15; p = 4, N = K = 3..10; по 10 пар;
пара k строится от зерна 7000003*N + 104729*p + k.

    python exp_translation.py [файл результатов] [-j число процессов]
"""
import argparse
import json
import os
import random
import time
from multiprocessing import Pool

from fsmnet import Net, rand_fsm, mb, length, second_order_sample, third_order_sample, translate
from wmethod import w_test


def one(args):
    N, p, k = args
    seed = 7000003 * N + 104729 * p + k
    rng = random.Random(seed)
    A = rand_fsm(N, p, p, rng); B = rand_fsm(N, p, p, rng)
    net = Net(A, B, p)
    RT = [net.mutant(f) for f in net.first_order(('T',))]
    m1 = [M for M in RT if net.distinguishable(M)]
    m2 = [M for f, M in second_order_sample(net, rng, 1000)]
    m3 = third_order_sample(net, rng, 1000)
    T1, _ = mb(net, net.cpp_tour(), RT, ext=False)
    # B в изоляции: головная компонента - тождественный автомат с одним состоянием
    ident = ([[0] * p], [list(range(p))])
    nb = Net(ident, B, p)
    if N <= 4:
        TB = w_test(nb, N); kind = 'W'
    else:
        RTb = [nb.mutant(f) for f in nb.first_order(('T',))]
        TB1, _ = mb(nb, nb.cpp_tour(), RTb, ext=False)
        so = [M for f, M in second_order_sample(nb, rng, 1000)]
        TB, _ = mb(nb, TB1, so, ext=False); kind = 'MB12'
    TT, full, sym = [], 0, 0
    for beta in TB:
        alpha = translate(A, beta, p)
        sym += len(alpha); full += (len(alpha) == len(beta))
        if alpha: TT.append(alpha)
    rec = dict(N=N, p=p, k=k, seed=seed, kind=kind, r=net.r, iso_len=length(TB), iso_cnt=len(TB), tt_len=length(TT),
               seq=full / len(TB) if TB else 1.0, sym=sym / max(1, sum(len(b) for b in TB)), MB1_len=length(T1))
    for lvl, ms in (('1', m1), ('2', m2), ('3', m3)):
        rec['n' + lvl] = len(ms)
        rec['TT' + lvl] = sum(net.detects(TT, M) for M in ms)
        rec['MB1_' + lvl] = sum(net.detects(T1, M) for M in ms)
    return rec


def main():
    ap = argparse.ArgumentParser(description='Трансляция теста изолированной компоненты: 350 пар случайных автоматов')
    ap.add_argument('out', nargs='?', default=os.path.join('results', 'translation_350.json'), help='файл результатов (по умолчанию %(default)s)')
    ap.add_argument('-j', '--procs', type=int, default=os.cpu_count() or 1, help='число процессов (по умолчанию - число ядер)')
    a = ap.parse_args()
    out, procs = a.out, a.procs
    jobs = [(N, p, k) for p, Ns in ((2, range(2, 16)), (3, range(3, 16)), (4, range(3, 11)))
            for N in Ns for k in range(10)]
    t0 = time.time()
    with Pool(procs) as pool:
        recs = pool.map(one, jobs, chunksize=1)
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump(recs, f, indent=0)
    print(f'пар: {len(recs)}, время {time.time() - t0:.0f} с, результаты: {out}')


if __name__ == '__main__':
    main()
