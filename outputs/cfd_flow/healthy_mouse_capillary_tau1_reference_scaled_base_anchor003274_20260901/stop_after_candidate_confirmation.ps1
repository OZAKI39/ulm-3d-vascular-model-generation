$ErrorActionPreference = "Stop"
$runRoot = "E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_flow\healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901"
$statePath = Join-Path $runRoot "qc\reference_scaled_base_run_state.json"
$runtime = "\\wsl.localhost\Ubuntu\home\lzy\u3da\tau1_reference_scaled_base_20260901"
$restart = Join-Path $runtime "restart"
$log = Join-Path $runRoot "confirmation_stop_monitor.log"
$shortWindow = 119751

while ($true) {
    if (-not (Test-Path -LiteralPath $statePath)) {
        Start-Sleep -Seconds 5
        continue
    }
    try {
        $state = Get-Content -LiteralPath $statePath -Raw | ConvertFrom-Json
    }
    catch {
        Start-Sleep -Seconds 2
        continue
    }
    if ($state.status -ne "IN_PROGRESS") {
        Add-Content -LiteralPath $log -Value "monitor exit: state=$($state.status)"
        exit 0
    }
    if ($null -eq $state.candidate_iteration) {
        Start-Sleep -Seconds 5
        continue
    }
    $candidate = [int64]$state.candidate_iteration
    $target = $candidate + $shortWindow
    $binary = Join-Path $restart "tau1_reference_scaled_base_$target.lsb"
    $header = Join-Path $restart "tau1_reference_scaled_base_header_$target.lua"
    if ((Test-Path -LiteralPath $binary) -and (Test-Path -LiteralPath $header)) {
        $firstSize = (Get-Item -LiteralPath $binary).Length
        Start-Sleep -Seconds 2
        $secondSize = (Get-Item -LiteralPath $binary).Length
        if ($firstSize -eq 27712640 -and $secondSize -eq $firstSize) {
            New-Item -ItemType File -Path (Join-Path $runtime "stop") -Force | Out-Null
            Add-Content -LiteralPath $log -Value "confirmation checkpoint $target complete; stop requested"
            exit 0
        }
    }
    Start-Sleep -Seconds 5
}
