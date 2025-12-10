$datetime = Get-Date -Format 'yyyy-MM-dd HH:mm:ss'
$out = New-Object System.Collections.Generic.List[string]
$out.Add('# Summary')
$out.Add('Run summarize.bat to update this file.')
$out.Add("## Last Updated: $datetime")
$out.Add('## Author: HuangZy')
$out.Add('[TOC]')

$sections = 'DataStructure','Geometry','Graph','Math'
foreach($sec in $sections){
  if(Test-Path $sec){
    Get-ChildItem -Recurse -Path $sec -Filter *.cpp | ForEach-Object {
      $out.Add('')
      $out.Add("## $($_.Name)")
      $out.Add('```cpp')
      $out.Add((Get-Content -Raw -Encoding UTF8 $_.FullName))
      $out.Add('```')
    }
  }
}

$out | Out-File -FilePath README.md -Encoding utf8

# powershell -ExecutionPolicy Bypass -File .\generate-readme.ps1