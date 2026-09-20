$h = Get-Content 'ar5iv_2110.12206.html' -Raw -Encoding UTF8
$pat = 'mathtable\b|marith[^o]|a11bits|a12bits|zltid|ltone[^s]|ltzero|l\(.?\)?M|half.?t0|mathname|\\alpha\b'
$ms = [regex]::Matches($h, $pat, 'IgnoreCase')
Write-Output ('MATCH COUNT=' + $ms.Count)
foreach($m in $ms){
  $s = [Math]::Max(0, $m.Index-140)
  $L = [Math]::Min(330, $h.Length-$s)
  $ctx = $h.Substring($s,$L)
  $ctx = $ctx -replace '\s+',' '
  Write-Output ('@' + $m.Index + ' :: ' + $ctx)
}
