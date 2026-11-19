#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Сводка результатов экспериментов в виде чисел, приведенных в статье.
Читает файлы из каталога results (по умолчанию) или из каталога, указанного аргументом:
    python report.py [каталог]
"""
import argparse
import json
import os
from collections import defaultdict
from statistics import mean

R = 'results'


def load(name):
    with open(os.path.join(R, name), encoding='utf-8') as f:
        return json.load(f)


def pct(x, d=1):
    return f'{100 * x:.{d}f}'.replace('.', ',') + ' %'


def num(x):
    return f'{x:,}'.replace(',', ' ')


def random_series():
    recs = load('random_436.json')
    print('=== Дополнительная серия: случайные автоматы (п. 3.1) ===')
    g = defaultdict(int)
    for r in recs: g[(r['p'], r['N'])] += 1
    print('пар:', len(recs), ' по (p, N=K):', ', '.join(f'({p}, {n}): {c}' for (p, n), c in sorted(g.items())))
    full = sum(r['MBC_M1'] == 1.0 for r in recs)
    print(f'MBC обнаруживает все обнаружимые мутанты 1-го порядка (переходы и выходы): в {full} из {len(recs)} пар')
    full = sum(r['MB1_M1'] == 1.0 for r in recs)
    print(f'MB1 обнаруживает все обнаружимые мутанты 1-го порядка: в {full} из {len(recs)} пар')
    n2 = sum(r['n2'] for r in recs)
    for name in ('MBC', 'MB1'):
        det = sum(r[name + '_2nd'] * r['n2'] for r in recs)
        print(f'{name}: обнаружено {pct(det / n2, 2)} из {num(round(n2))} обнаружимых мутантов 2-го порядка (выборки)')
    print('длина MBC относительно MB1, в среднем по парам:', pct(mean(r['MBC'] / r['MB1'] for r in recs)))
    by = defaultdict(list)
    for r in recs: by[(r['p'], r['N'])].append(r['MBC'] / r['MB1'])
    for key in sorted(by):
        if key[1] >= 15:
            print(f'  p = {key[0]}, N = K = {key[1]}: {pct(mean(by[key]))}')
    print()


def translation_series():
    recs = load('translation_350.json')
    print('=== Трансляция теста изолированной компоненты: случайные автоматы (п. 3.2) ===')
    print('пар:', len(recs))
    print('доля последовательностей, транслируемых целиком:', pct(mean(r['seq'] for r in recs)))
    print('доля входных символов, доходящих до B:', pct(mean(r['sym'] for r in recs)))

    def cov(rs, pref):
        return [sum(r[pref + l] for r in rs) / sum(r['n' + l] for r in rs) for l in '123']
    print('полнота транслированного теста (1-й, 2-й, 3-й порядок):', ', '.join(pct(x) for x in cov(recs, 'TT')))
    p2 = [r for r in recs if r['p'] == 2]
    print('  то же при p = 2:', ', '.join(pct(x) for x in cov(p2, 'TT')))
    print('полнота MB1 на тех же мутантах:', ', '.join(pct(x, 2) for x in cov(recs, 'MB1_')))
    full = sum(r['TT1'] == r['n1'] for r in recs)
    print(f'транслированный тест обнаружил все мутанты 1-го порядка в {full} из {len(recs)} пар')
    print()


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Сводка результатов экспериментов')
    ap.add_argument('dir', nargs='?', default='results', help='каталог с результатами (по умолчанию %(default)s)')
    R = ap.parse_args().dir
    random_series()
    translation_series()
