param([string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$projectDirectory = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($OutputDirectory)) {
    $OutputDirectory = Join-Path $env:TEMP 'hmb-explanation-20261104-build'
}
$buildDirectory = [System.IO.Path]::GetFullPath($OutputDirectory)
if ($buildDirectory -eq [System.IO.Path]::GetFullPath($projectDirectory)) {
    throw 'Use a separate build directory for auxiliary files.'
}
New-Item -ItemType Directory -Force -Path $buildDirectory | Out-Null
Push-Location $projectDirectory
try {
    & xelatex -recorder -synctex=1 -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$buildDirectory" hmb-explanation.tex > (Join-Path $buildDirectory 'pass1.txt') 2>&1
    if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass 1 failed; see $buildDirectory\pass1.txt" }
    Copy-Item -LiteralPath (Join-Path $projectDirectory 'citelist.bib') -Destination $buildDirectory -Force
    Copy-Item -LiteralPath (Join-Path $projectDirectory 'gbt7714-unsrt.bst') -Destination $buildDirectory -Force
    Push-Location $buildDirectory
    try {
        & bibtex hmb-explanation > (Join-Path $buildDirectory 'bibtex.txt') 2>&1
        if ($LASTEXITCODE -ne 0) { throw "BibTeX failed; see $buildDirectory\bibtex.txt" }
    }
    finally { Pop-Location }
    foreach ($pass in 2,3) {
        & xelatex -recorder -synctex=1 -interaction=nonstopmode -halt-on-error -file-line-error "-output-directory=$buildDirectory" hmb-explanation.tex > (Join-Path $buildDirectory "pass$pass.txt") 2>&1
        if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass $pass failed; see $buildDirectory\pass$pass.txt" }
    }
    Copy-Item -LiteralPath (Join-Path $buildDirectory 'hmb-explanation.pdf') -Destination (Join-Path $projectDirectory 'hmb-explanation.pdf') -Force
    Copy-Item -LiteralPath (Join-Path $buildDirectory 'hmb-explanation.synctex.gz') -Destination (Join-Path $projectDirectory 'hmb-explanation.synctex.gz') -Force
    Write-Output "BUILD SUCCESS: $projectDirectory\hmb-explanation.pdf"
    Write-Output "Build directory: $buildDirectory"
}
finally { Pop-Location }
