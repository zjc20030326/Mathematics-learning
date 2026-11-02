param([string]$OutputDirectory = '')
$ErrorActionPreference = 'Stop'
$projectDirectory = $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($OutputDirectory)) {
    $OutputDirectory = Join-Path $env:TEMP 'hmb-explanation-build'
}
$buildDirectory = [System.IO.Path]::GetFullPath($OutputDirectory)
New-Item -ItemType Directory -Force -Path $buildDirectory | Out-Null
Copy-Item -LiteralPath (Join-Path $projectDirectory 'citelist.bib') -Destination $buildDirectory -Force

Push-Location $projectDirectory
try {
    & latexmk -xelatex -interaction=nonstopmode -halt-on-error "-output-directory=$buildDirectory" hmb-explanation.tex
    if ($LASTEXITCODE -ne 0) { throw "LaTeXmk compilation failed; check output." }
    Copy-Item -LiteralPath (Join-Path $buildDirectory 'hmb-explanation.pdf') -Destination (Join-Path $projectDirectory 'hmb-explanation.pdf') -Force
    Write-Output "BUILD SUCCESS: $projectDirectory\hmb-explanation.pdf"
    Write-Output "Build directory: $buildDirectory"
}
finally { Pop-Location }
