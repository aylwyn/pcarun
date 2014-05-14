#!/usr/bin/awk -f

BEGIN{OFS = "\t"}

!/^#/ && $3=="."{$3 = $1":"$2}

{print}
