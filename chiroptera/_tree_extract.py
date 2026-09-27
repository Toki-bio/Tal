import re
import sys

h = open(sys.argv[1], encoding="utf-8").read()
i = h.find('<svg viewBox="0 0 1280')
svg = h[i:h.find("</svg>", i) + 6]
svg = re.sub(r'width="100%" style="[^"]*"', 'width="1280"', svg, count=1)
open(sys.argv[2], "w", encoding="utf-8").write(svg)
