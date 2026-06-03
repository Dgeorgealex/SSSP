#!/usr/bin/env python3
"""
Create an EMD-style instance from two PGM images.

Usage:
  ./create_emd_distance.py imgA.pgm imgB.pgm out.txt [--8]

Output format (compatible with src/cycle_cancelling.cpp::read_graph):
  n m
  b0 b1 b2 ... b_{n-1}
  u v
  u v
  ... (m lines of directed edges)

Each node is a pixel (row-major id = y*width + x).
Edges: 4-neighbour grid by default (add both directions). Use --8 to add diagonal neighbours as well.
Balance at a node = imgA - imgB (signed integer). The script will print a warning if total balance != 0.
"""

import sys
from pathlib import Path
import argparse


def read_pgm(path):
    with open(path, 'rb') as f:
        header = f.readline().strip()
        if header not in (b'P5', b'P2'):
            raise ValueError(f'Unsupported PGM format: {header!r}')
        # skip comments
        def read_token():
            while True:
                line = f.readline()
                if not line:
                    return b''
                line = line.strip()
                if line.startswith(b'#') or len(line) == 0:
                    continue
                for tok in line.split():
                    yield tok
        tokgen = read_token()
        tokens = []
        # read width, height, maxval
        try:
            w = int(next(tokgen))
            h = int(next(tokgen))
            maxv = int(next(tokgen))
        except StopIteration:
            raise ValueError('Invalid PGM header')
        # if P5, read binary pixels; if P2 we need to continue token parsing
        pixels = []
        if header == b'P5':
            # consume one byte if the next byte is newline
            # f is already at the next byte after the last token; ensure we're at pixel start
            # read remaining whitespace until a single byte remains
            # In practice after reading tokens with token generator we can't easily resume.
            # So reparse header manually for robustness.
            f.seek(0)
            content = f.read()
            # simple parser: split header parts
            parts = content.split()
            # parts: [P5, w, h, maxv, <binary pixels...>]
            header_len = 4
            raw = b' '.join(parts[header_len:])
            # For P5 with maxv < 256 each pixel is one byte; else two bytes big-endian
            if maxv < 256:
                expected = w * h
                raw_pixels = raw[:expected]
                if len(raw_pixels) < expected:
                    raise ValueError('Not enough pixel data')
                pixels = list(raw_pixels)
            else:
                # 2 bytes per sample (big endian)
                expected = w * h * 2
                raw_pixels = raw[:expected]
                if len(raw_pixels) < expected:
                    raise ValueError('Not enough pixel data')
                pixels = [ (raw_pixels[i]<<8) | raw_pixels[i+1] for i in range(0, expected, 2) ]
        else:
            # P2 ASCII; reuse token generator by reopening as text
            with open(path, 'rt') as ft:
                toks = []
                for line in ft:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    toks.extend(line.split())
                # toks[0] = P2, toks[1]=w, toks[2]=h, toks[3]=maxv, rest pixels
                if len(toks) < 4 + w*h:
                    raise ValueError('Not enough pixel tokens')
                pixels = [int(x) for x in toks[4:4+w*h]]
        if len(pixels) != w*h:
            raise ValueError(f'Pixel count mismatch: got {len(pixels)}, expected {w*h}')
        return w, h, maxv, pixels


def build_edges(w, h, use8=False):
    edges = []
    for y in range(h):
        for x in range(w):
            u = y*w + x
            # 4-neighbour
            if x+1 < w:
                v = y*w + (x+1)
                edges.append((u, v))
                edges.append((v, u))
            if y+1 < h:
                v = (y+1)*w + x
                edges.append((u, v))
                edges.append((v, u))
            if use8:
                if x+1 < w and y+1 < h:
                    v = (y+1)*w + (x+1)
                    edges.append((u, v)); edges.append((v, u))
                if x-1 >= 0 and y+1 < h:
                    v = (y+1)*w + (x-1)
                    edges.append((u, v)); edges.append((v, u))
    return edges


def main():
    p = argparse.ArgumentParser()
    p.add_argument('imgA', type=Path)
    p.add_argument('imgB', type=Path)
    p.add_argument('out', type=Path, nargs='?', default=Path('graph.txt'))
    p.add_argument('--8', dest='use8', action='store_true', help='use 8-neighbour connectivity')
    args = p.parse_args()

    if not args.imgA.exists() or not args.imgB.exists():
        print('Input image(s) not found')
        sys.exit(1)

    w1,h1,max1,pix1 = read_pgm(str(args.imgA))
    w2,h2,max2,pix2 = read_pgm(str(args.imgB))
    if (w1,h1) != (w2,h2):
        print('Image sizes differ')
        sys.exit(1)

    w,h = w1,h1
    n = w*h
    # balances: imgA - imgB
    balances = [int(a) - int(b) for a,b in zip(pix1, pix2)]
    s = sum(balances)
    if s != 0:
        print(f'Warning: total balance != 0 ({s}). EMD requires total supply == total demand.')

    edges = build_edges(w,h,use8=args.use8)
    m = len(edges)

    # write output
    with open(args.out, 'w') as out:
        out.write(f"{n} {m}\n")
        for b in balances:
            out.write(f"{b}\n")
        for (u,v) in edges:
            out.write(f"{u} {v}\n")

    print(f'Wrote {args.out} with n={n}, m={m}')

if __name__ == '__main__':
    main()
