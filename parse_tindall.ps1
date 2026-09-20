Set-Location 'D:\MUBs in 6-dimension'

function Join-Array($items) {
  if ($null -eq $items) { return '' }
  if ($items -is [array]) { return ($items | ForEach-Object { [string]$_ }) -join ';' }
  return [string]$items
}

# OpenAlex author works
$oa = Get-Content -Raw 'tindall_openalex.json' | ConvertFrom-Json
$out = @()
foreach ($w in $oa.results) {
  $venue = $null
  if ($null -ne $w.primary_location.source) { $venue = $w.primary_location.source.display_name }
  $authors = @()
  foreach ($a in $w.authorships) {
    $aid = ''
    if ($null -ne $a.author.id) { $aid = $a.author.id }
    $authors += ($a.author.display_name + ' [' + $aid + ']')
  }
  $ids = @()
  foreach ($k in $w.ids.PSObject.Properties) { $ids += $k.Name + '=' + $k.Value }
  $out += [pscustomobject]@{
    openalex_id = $w.id
    doi = $w.doi
    title = $w.title
    publication_date = $w.publication_date
    publication_year = $w.publication_year
    type = $w.type
    venue = $venue
    authors = $authors
    ids = $ids
    landing = $w.primary_location.landing_page_url
    oa_url = $w.open_access.oa_url
    cited = $w.cited_by_count
  }
}
$out | Export-Csv 'tindall_openalex_records.csv' -NoTypeInformation -Encoding utf8

# Crossref query response: filter exact Joseph Tindall authorship and retain all records
$cr = Get-Content -Raw 'tindall_crossref.json' | ConvertFrom-Json
$crout = @()
foreach ($w in $cr.message.items) {
  $has = $false
  foreach ($a in $w.author) {
    $name = (($a.given + ' ' + $a.family).Trim())
    if ($name -eq 'Joseph Tindall' -or $name -eq 'Tindall, Joseph') { $has = $true }
  }
  if (-not $has) { continue }
  $authors = @()
  foreach ($a in $w.author) { $authors += (($a.given + ' ' + $a.family).Trim()) }
  $rels = @()
  if ($null -ne $w.relation) { foreach ($r in $w.relation.PSObject.Properties) { $rels += $r.Name + '=' + ($r.Value | ConvertTo-Json -Compress) } }
  $crout += [pscustomobject]@{
    doi = $w.DOI
    title = Join-Array $w.title
    published = Join-Array $w.published
    created = Join-Array $w.created
    issued = Join-Array $w.issued
    type = $w.type
    publisher = Join-Array $w.publisher
    container = Join-Array $w['container-title']
    volume = Join-Array $w.volume
    issue = Join-Array $w.issue
    page = Join-Array $w.page
    article_number = Join-Array $w['article-number']
    authors = $authors
    url = $w.URL
    relation = $rels
  }
}
$crout | Export-Csv 'tindall_crossref_records.csv' -NoTypeInformation -Encoding utf8

# ORCID JSON record
$oc = Get-Content -Raw 'tindall_orcid.json' | ConvertFrom-Json
$ocout = @()
foreach ($g in $oc.'orcid-record'.activities.group) {
  foreach ($w in $g.works) {
    $dois = @()
    foreach ($eid in $w.'external-ids'.'external-id') {
      if ($eid.'external-id-type' -eq 'doi') { $dois += $eid.'external-id-value' }
    }
    $auth = @()
    foreach ($aw in $w.'author-credit') { $auth += ($aw.'given-names' + ' ' + $aw.'family-name') }
    $ext = @()
    foreach ($eid in $w.'external-ids'.'external-id') { $ext += $eid.'external-id-type' + '=' + $eid.'external-id-value' }
    $ocout += [pscustomobject]@{
      path = $w.path
      title = $w.title.title
      type = $w.type
      doi = $dois
      external_ids = $ext
      publication_date = ($w.'publication-date'.year + '-' + $w.'publication-date'.month + '-' + $w.'publication-date'.day)
      authors = $auth
      source = $w.'source'.'source-name'
    }
  }
}
$ocout | Export-Csv 'tindall_orcid_records.csv' -NoTypeInformation -Encoding utf8

# arXiv XML
[xml]$ax = Get-Content -Raw 'tindall_arxiv.xml'
$ns = @{ a='http://www.w3.org/2005/Atom'; x='http://arxiv.org/schemas/atom' }
$ao = @()
foreach ($n in $ax.feed.entry) {
  $authors = @()
  foreach ($au in $n.author) { $authors += $au.name }
  $cats = @()
  foreach ($c in $n.category) { $cats += $c.'#text' }
  $dois = @()
  foreach ($d in $n.'arxiv:doi') { $dois += $d.'#text' }
  $refs = @()
  foreach ($r in $n.'arxiv:journal_ref') { $refs += $r.'#text' }
  $comments = @()
  foreach ($cmt in $n.'arxiv:comment') { $comments += $cmt.'#text' }
  $ao += [pscustomobject]@{
    id = $n.id
    title = $n.title
    updated = $n.updated
    published = $n.published
    authors = $authors
    journal_ref = $refs
    doi = $dois
    comment = $comments
    categories = $cats
    link = $n.link.href
  }
}
$ao | Export-Csv 'tindall_arxiv_records.csv' -NoTypeInformation -Encoding utf8

"OpenAlex $($oa.meta.count)"
"Crossref matching $($crout.Count)"
"ORCID $($ocout.Count)"
"arXiv $($ao.Count)"
