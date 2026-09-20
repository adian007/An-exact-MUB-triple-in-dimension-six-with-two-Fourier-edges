$ProgressPreference='SilentlyContinue'
try {
  $r = Invoke-WebRequest -Uri 'https://ar5iv.labs.arxiv.org/html/2110.12206' -UseBasicParsing -TimeoutSec 30
  Write-Output ('STATUS=' + $r.StatusCode)
  Write-Output ('LEN=' + $r.Content.Length)
  $r.Content | Out-File -Encoding utf8 'D:\MUBs in 6-dimension\_spawn_cache\paper_2110.12206.html'
  Write-Output 'SAVED'
} catch {
  Write-Output ('ERR: ' + $_.Exception.Message)
}
