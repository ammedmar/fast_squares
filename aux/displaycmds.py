'''
Reads defined commands in and creates tex file to render them all

'''

import re
import os

skip = {'pdfEinfty', 'anibal'}

with open('usualcmds.tex', 'r') as f:
    contents = f.read()

    pattern1 = re.compile(r'[dr]{(.+?)}{')
    matches1 = pattern1.finditer(contents)

    pattern2 = re.compile(r'[dr]{(.+?)}\[(\d)\]')
    matches2 = pattern2.finditer(contents)

    with open('displaycmds.tex', 'w') as fp:
        for m in matches1:
            if all([w not in m[1] for w in skip]):
                fp.write(f'\\verb|{m[1]}|, ${m[1]}$ ; \n')
        for m in matches2:
            if m[2] == '1':
                fp.write(f'\\verb|{m[1]}{{x}}|, ${m[1]}{{x}}$ ; \n')
            if m[2] == '2':
                fp.write(f'\\verb|{m[1]}{{x}}{{y}}|, ${m[1]}{{x}}{{y}}$ ; \n')
