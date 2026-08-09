[CmdletBinding()]
param(
    [string]$Executable = (Join-Path $PSScriptRoot '..\build\Bin\Release\GREAT_PVT.exe'),
    [string]$ConfigDirectory = (Join-Path $PSScriptRoot '..\build\parallel30'),
    [string[]]$Config = @(
        'UPD_NONE.xml',
        'UPD_PARAMETER.xml',
        'UPD_CONSTRAINT.xml',
        'UPD_CONSTRAINT_GODN.xml',
        'OSB_NONE.xml',
        'OSB_PARAMETER.xml',
        'OSB_CONSTRAINT.xml',
        'OSB_CONSTRAINT_HARB.xml'
    ),
    [int]$MaxParallel = 6,
    [switch]$DryRun
)

# Run independent FGO processes in parallel; RAW Ceres workers remain serialized safely.
$repoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$exePath = (Resolve-Path $Executable).Path
$configRoot = (Resolve-Path $ConfigDirectory).Path
if ($MaxParallel -lt 1) { throw 'MaxParallel must be at least 1.' }

$jobs = foreach ($configName in $Config) {
    $configPath = (Resolve-Path (Join-Path $configRoot $configName)).Path
    [xml]$xml = Get-Content -LiteralPath $configPath -Raw
    $begin = [DateTime]::Parse($xml.config.gen.beg, [Globalization.CultureInfo]::InvariantCulture)
    $end = [DateTime]::Parse($xml.config.gen.end, [Globalization.CultureInfo]::InvariantCulture)
    $duration = $end - $begin
    if ($duration.TotalMinutes -gt 30.0001) {
        throw "$configName covers $($duration.TotalMinutes.ToString('0.###')) minutes; expected at most 30 minutes."
    }

    # Config paths and input files are relative to their matching sample-data directory.
    $isOsb = $xml.config.inputs.bias -match 'OSB|WUM0MGXRAP'
    $dataDirName = if ($isOsb) { 'PPPFLT_2023305_OSB' } else { 'PPPFLT_2023305' }
    $workDir = (Resolve-Path (Join-Path $repoRoot "sample_data\$dataDirName")).Path
    $relativeConfig = [IO.Path]::GetRelativePath($workDir, $configPath)
    $logStem = [IO.Path]::GetFileNameWithoutExtension($configName)
    [PSCustomObject]@{
        Name = $logStem
        Config = $configName
        ConfigArgument = $relativeConfig
        WorkDir = $workDir
        Stdout = Join-Path $configRoot "$logStem.parallel.stdout.log"
        Stderr = Join-Path $configRoot "$logStem.parallel.stderr.log"
    }
}

if ($DryRun) {
    $jobs | Format-Table Name, Config, ConfigArgument, WorkDir -AutoSize
    exit 0
}

$pending = [Collections.Generic.Queue[object]]::new()
$jobs | ForEach-Object { $pending.Enqueue($_) }
$active = @{}
$results = [Collections.Generic.List[object]]::new()

while ($pending.Count -gt 0 -or $active.Count -gt 0) {
    while ($pending.Count -gt 0 -and $active.Count -lt $MaxParallel) {
        $job = $pending.Dequeue()
        $process = Start-Process -FilePath $exePath `
            -ArgumentList @('-x', $job.ConfigArgument) `
            -WorkingDirectory $job.WorkDir `
            -RedirectStandardOutput $job.Stdout `
            -RedirectStandardError $job.Stderr `
            -PassThru
        $active[$process.Id] = [PSCustomObject]@{ Process = $process; Job = $job; Started = Get-Date }
        Write-Host ("START {0} (pid {1})" -f $job.Name, $process.Id)
    }

    Start-Sleep -Milliseconds 500
    foreach ($id in @($active.Keys)) {
        $entry = $active[$id]
        if (-not $entry.Process.HasExited) { continue }
        $entry.Process.WaitForExit()
        $elapsed = ((Get-Date) - $entry.Started).TotalSeconds
        $results.Add([PSCustomObject]@{
            Name = $entry.Job.Name
            ExitCode = $entry.Process.ExitCode
            WallSeconds = [Math]::Round($elapsed, 3)
            Stdout = $entry.Job.Stdout
            Stderr = $entry.Job.Stderr
        })
        $active.Remove($id)
        Write-Host ("DONE  {0} exit={1} wall={2}s" -f $entry.Job.Name, $entry.Process.ExitCode, [Math]::Round($elapsed, 3))
    }
}

$results | Sort-Object Name | Format-Table Name, ExitCode, WallSeconds -AutoSize
$failed = @($results | Where-Object ExitCode -ne 0)
if ($failed.Count -gt 0) {
    throw ("{0} FGO process(es) failed: {1}" -f $failed.Count, (($failed | ForEach-Object Name) -join ', '))
}
