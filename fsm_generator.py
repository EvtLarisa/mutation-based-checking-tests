#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Генератор случайных пар конечных автоматов A и B для экспериментов
с последовательной композицией A∙B (статья «Исследование мутационных
проверяющих тестов для компонентов автоматных сетей»).

Каждый автомат детерминированный, полностью определенный, инициальный
(начальное состояние 0). Функции переходов и выходов выбираются
равномерно и независимо для каждой пары (состояние, входной символ).

Композиция: внешние входы I подаются на A, выход A (внутренний канал U)
является входом B, выход B (канал O) наблюдаем; |I| = |U| = |O| = p.

Генерация воспроизводима: пара с номером k для параметров (N, K, p)
строится от зерна seed = 1000003*N + 7919*p + k, сначала A, потом B,
от одного и того же генератора случайных чисел. С этими зернами были
получены автоматы дополнительной серии экспериментов (ключ --article,
436 пар). Пары серии с трансляцией теста (ключ --translation, 350 пар)
строятся по тем же правилам от зерна 7000003*N + 104729*p + k.

Только стандартная библиотека Python (3.8+).

Примеры:
  python fsm_generator.py --N 5 --K 5 --p 2 --pairs 10 --out pairs
  python fsm_generator.py --article --out article_pairs
  python fsm_generator.py --translation --out translation_pairs
  python fsm_generator.py --N 3 --K 3 --p 2 --pairs 1 --print
"""
import argparse
import json
import os
import random
from collections import deque


def rand_fsm(n, nin, nout, rng):
    """Случайный автомат с n состояниями, nin входами и nout выходами.
    delta[s][i] - следующее состояние, lam[s][i] - выходной символ."""
    delta = [[rng.randrange(n) for _ in range(nin)] for _ in range(n)]
    lam = [[rng.randrange(nout) for _ in range(nin)] for _ in range(n)]
    return delta, lam


def pair_seed(N, p, k, scheme='article'):
    if scheme == 'translation':
        return 7000003 * N + 104729 * p + k
    return 1000003 * N + 7919 * p + k


def generate_pair(N, K, p, seed):
    """Пара (A, B): A имеет N состояний, B - K состояний, все алфавиты мощности p."""
    rng = random.Random(seed)
    A = rand_fsm(N, p, p, rng)
    B = rand_fsm(K, p, p, rng)
    return A, B


def composition(A, B, p):
    """Достижимая часть автомата-композиции A∙B (начальное состояние (0, 0)).
    Возвращает список состояний (a, b) и таблицы переходов/выходов по входам I."""
    (dA, lA), (dB, lB) = A, B
    idx = {(0, 0): 0}
    states = [(0, 0)]
    q = deque([(0, 0)])
    delta, lam = [], []
    while q:
        a, b = q.popleft()
        row_d, row_l = [], []
        for i in range(p):
            u = lA[a][i]
            t = (dA[a][i], dB[b][u])
            if t not in idx:
                idx[t] = len(states)
                states.append(t)
                q.append(t)
            row_d.append(idx[t])
            row_l.append(lB[b][u])
        delta.append(row_d)
        lam.append(row_l)
    return states, delta, lam


def fsm_dict(fsm, inputs, outputs):
    delta, lam = fsm
    return {'states': len(delta), 'inputs': inputs, 'outputs': outputs, 'initial': 0,
            'delta': delta, 'lambda': lam}


def fsm_text(name, fsm, in_name, out_name):
    """Текстовая таблица: строка «состояние вход следующее_состояние выход»."""
    delta, lam = fsm
    n, p = len(delta), len(delta[0])
    lines = [f'# {name}: состояний {n}, входной алфавит {in_name} = {{0..{p - 1}}}, '
             f'выходной алфавит {out_name} = {{0..{p - 1}}}, начальное состояние 0',
             '# s i s\' o']
    for s in range(n):
        for i in range(p):
            lines.append(f'{s} {i} {delta[s][i]} {lam[s][i]}')
    return '\n'.join(lines)


def article_grid():
    """Параметры дополнительной серии экспериментов: (N, K, p, число пар) - всего 436 пар."""
    grid = []
    for p in (2, 3, 4):
        for n in range(3, 16, 2):
            grid.append((n, n, p, 20))
    grid += [(20, 20, 2, 5), (30, 30, 2, 5), (40, 40, 2, 3), (50, 50, 2, 3)]
    return grid


def translation_grid():
    """Параметры серии с трансляцией теста: (N, K, p, число пар) - всего 350 пар."""
    return ([(n, n, 2, 10) for n in range(2, 16)] + [(n, n, 3, 10) for n in range(3, 16)] +
            [(n, n, 4, 10) for n in range(3, 11)])


def save_pair(out, N, K, p, k, seed, A, B):
    states, _, _ = composition(A, B, p)
    base = os.path.join(out, f'N{N}_K{K}_p{p}_{k:02d}')
    rec = {'N': N, 'K': K, 'p': p, 'pair': k, 'seed': seed,
           'composition_states': len(states),
           'A': fsm_dict(A, 'I', 'U'), 'B': fsm_dict(B, 'U', 'O')}
    with open(base + '.json', 'w', encoding='utf-8') as f:
        json.dump(rec, f, ensure_ascii=False)
    with open(base + '.txt', 'w', encoding='utf-8') as f:
        f.write(f'# N = {N}, K = {K}, p = {p}, пара {k}, seed = {seed}, '
                f'достижимых состояний A∙B: {len(states)}\n')
        f.write(fsm_text('A', A, 'I', 'U') + '\n\n' + fsm_text('B', B, 'U', 'O') + '\n')
    return len(states)


def main():
    ap = argparse.ArgumentParser(description='Генератор случайных пар конечных автоматов для композиции A∙B')
    ap.add_argument('--N', type=int, help='число состояний автомата A')
    ap.add_argument('--K', type=int, help='число состояний автомата B (по умолчанию K = N)')
    ap.add_argument('--p', type=int, default=2, help='мощность алфавитов |I| = |U| = |O| (по умолчанию 2)')
    ap.add_argument('--pairs', type=int, default=10, help='число пар (по умолчанию 10)')
    ap.add_argument('--first', type=int, default=0, help='номер первой пары k (по умолчанию 0)')
    ap.add_argument('--article', action='store_true', help='все 436 пар дополнительной серии экспериментов')
    ap.add_argument('--translation', action='store_true', help='все 350 пар серии с трансляцией теста')
    ap.add_argument('--out', default='pairs', help='каталог для файлов .json и .txt')
    ap.add_argument('--print', action='store_true', help='только вывести автоматы на экран')
    a = ap.parse_args()

    scheme = 'article'
    if a.article:
        grid = article_grid()
    elif a.translation:
        grid = translation_grid(); scheme = 'translation'
    else:
        if a.N is None:
            ap.error('укажите --N (или --article, --translation)')
        grid = [(a.N, a.K if a.K is not None else a.N, a.p, a.pairs)]

    if not a.print:
        os.makedirs(a.out, exist_ok=True)
    total = 0
    for N, K, p, npairs in grid:
        sizes = []
        for k in range(a.first, a.first + npairs):
            seed = pair_seed(N, p, k, scheme)
            A, B = generate_pair(N, K, p, seed)
            if a.print:
                states, _, _ = composition(A, B, p)
                print(f'# N = {N}, K = {K}, p = {p}, пара {k}, seed = {seed}, достижимых состояний A∙B: {len(states)}')
                print(fsm_text('A', A, 'I', 'U'))
                print()
                print(fsm_text('B', B, 'U', 'O'))
                print()
            else:
                sizes.append(save_pair(a.out, N, K, p, k, seed, A, B))
            total += 1
        if sizes:
            word = 'пары' if npairs % 10 in (2, 3, 4) and npairs % 100 not in (12, 13, 14) else 'пар'
            print(f'N = {N}, K = {K}, p = {p}: {npairs} {word}, '
                  f'в среднем {sum(sizes) / len(sizes):.1f} достижимых состояний композиции')
    if not a.print:
        print(f'всего пар: {total}, файлы в каталоге {a.out}')


if __name__ == '__main__':
    main()
