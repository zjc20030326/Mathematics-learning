param([string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$projectDirectory = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($OutputDirectory)) {
    $OutputDirectory = Join-Path $env:TEMP 'hmb-six-chapters-build'
}
$buildDirectory = [System.IO.Path]::GetFullPath($OutputDirectory)
New-Item -ItemType Directory -Force -Path $buildDirectory | Out-Null
Push-Location $projectDirectory
try {
    & xelatex -recorder -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$buildDirectory" fduthesis.tex > (Join-Path $buildDirectory 'pass1.txt') 2>&1
    if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass 1 failed; see $buildDirectory\pass1.txt" }
    Copy-Item -LiteralPath (Join-Path $projectDirectory 'citelist.bib') -Destination $buildDirectory -Force
    Get-ChildItem -LiteralPath $projectDirectory -Filter '*.bst' | Copy-Item -Destination $buildDirectory -Force
    Push-Location $buildDirectory
    try {
        & bibtex fduthesis > (Join-Path $buildDirectory 'bibtex.txt') 2>&1
        if ($LASTEXITCODE -ne 0) { throw "BibTeX failed; see $buildDirectory\bibtex.txt" }
    }
    finally { Pop-Location }
    foreach ($pass in 2,3) {
        & xelatex -recorder -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$buildDirectory" fduthesis.tex > (Join-Path $buildDirectory "pass$pass.txt") 2>&1
        if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass $pass failed; see $buildDirectory\pass$pass.txt" }
    }
    Copy-Item -LiteralPath (Join-Path $buildDirectory 'fduthesis.pdf') -Destination (Join-Path $projectDirectory 'fduthesis.pdf') -Force
    Write-Output "BUILD SUCCESS: $projectDirectory\fduthesis.pdf"
    Write-Output "Build directory: $buildDirectory"
}
finally { Pop-Location }
