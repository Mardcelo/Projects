## Fetch module from Tim Menzies' gist linked above
from a12 import *

## Create a labeled array
B_norm = ["baseline"]
## Append B values from baseline measurements
B_norm.extend(B)
## Likewise for tweak measurements
T_norm = ["tweak"]
T_norm.extend(T)
## Create consolidated list
C = [B_norm, T_norm]
for rx in a12s(C,rev=True,enough=0.71): print(rx)
