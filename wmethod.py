# -*- coding: utf-8 -*-
"""
W-тест, полный относительно реализаций с не более чем m состояниями, для приведенной формы
автомата Spec (класс Net из fsmnet; для компоненты B в изоляции головная компонента - тождественный
автомат с одним состоянием).

  W-тест: Q · I^{<=m-n+1} · W, Q - покрытие состояний, W - характеризующее множество.
Последовательности, являющиеся префиксами других, отбрасываются.
"""
import itertools
from collections import deque


def reduce_spec(net):
    """Приведенная форма достижимой части Spec: (n, nxt, out) по классам эквивалентности; класс s0 = 0."""
    r, p = net.r, net.p
    out = [[net.lB[b][u] for (b, u) in net.lab[s]] for s in range(r)]
    part = [tuple(out[s]) for s in range(r)]
    while True:
        ids = {}
        c = [ids.setdefault(v, len(ids)) for v in part]
        new = [(c[s],) + tuple(c[net.nxt[s][i]] for i in range(p)) for s in range(r)]
        if len(set(new)) == len(set(part)):
            cls = c; break
        part = new
    order = {}
    for s in range(r):
        order.setdefault(cls[s], len(order))
    n = len(order)
    rn = [None] * n; ro = [None] * n
    for s in range(r):
        k = order[cls[s]]
        if rn[k] is None:
            rn[k] = [order[cls[net.nxt[s][i]]] for i in range(p)]
            ro[k] = out[s][:]
    return n, rn, ro


def state_cover(n, nxt, p):
    path = {0: ()}
    q = deque([0])
    while q:
        s = q.popleft()
        for i in range(p):
            t = nxt[s][i]
            if t not in path:
                path[t] = path[s] + (i,); q.append(t)
    return [path[s] for s in range(n)]


def separating(n, nxt, out, p):
    """кратчайшая различающая последовательность для каждой пары состояний (s < t)"""
    sep = {}
    for s in range(n):
        for t in range(s + 1, n):
            seen = {(s, t): ()}
            q = deque([(s, t)])
            found = None
            while q and found is None:
                a, b = q.popleft()
                for i in range(p):
                    if out[a][i] != out[b][i]:
                        found = seen[(a, b)] + (i,); break
                    na, nb = nxt[a][i], nxt[b][i]
                    key = (min(na, nb), max(na, nb))
                    if na != nb and key not in seen:
                        seen[key] = seen[(a, b)] + (i,); q.append(key)
            sep[(s, t)] = found
    return sep


def maximal_len(seqs):
    """(длина со сбросами, число последовательностей, сами последовательности) после удаления префиксов"""
    ss = sorted(set(''.join(chr(48 + i) for i in s) for s in seqs if len(s) > 0))
    tot = 0; cnt = 0; keep = []
    for i, a in enumerate(ss):
        if i + 1 < len(ss) and ss[i + 1].startswith(a):
            continue
        tot += len(a) + 1; cnt += 1; keep.append(a)
    return tot, cnt, keep


def _decode(keep):
    return [[ord(ch) - 48 for ch in a] for a in keep]


def w_test(net, m):
    """только W-тест (для приведенной формы с одним состоянием W = {пустое слово})"""
    p = net.p
    n, nxt, out = reduce_spec(net)
    Q = state_cover(n, nxt, p)
    sep = separating(n, nxt, out, p)
    W = sorted(set(sep.values())) or [()]
    mids = [()]
    for L in range(1, m - n + 2):
        mids += list(itertools.product(range(p), repeat=L))
    Wt = [q + x + w for q in Q for x in mids for w in W]
    wl, wc, wk = maximal_len(Wt)
    return _decode(wk)
