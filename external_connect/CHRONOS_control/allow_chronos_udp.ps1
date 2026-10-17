# Run once as administrator. Scope: this camera, this Ethernet address, UDP 5000-5001.
$ErrorActionPreference = 'Stop'
$ruleName = 'Chronos14-192.168.0.10-UDP5000'
$logFile = Join-Path $PSScriptRoot 'firewall_setup.log'
try {
    $adapter = Get-NetIPAddress -AddressFamily IPv4 -IPAddress '192.168.0.101'
    if (-not $adapter) { throw 'The camera Ethernet address 192.168.0.101 is not present.' }
    $existing = Get-NetFirewallRule -Name $ruleName -ErrorAction SilentlyContinue
    if ($existing) {
        'Rule already exists. No changes made.' | Out-File -LiteralPath $logFile -Encoding utf8
    } else {
        New-NetFirewallRule -Name $ruleName -DisplayName 'Chronos 1.4 video from 192.168.0.10 only' `
            -Direction Inbound -Action Allow -Enabled True -Profile Public -Protocol UDP `
            -LocalPort 5000-5001 -LocalAddress '192.168.0.101' -RemoteAddress '192.168.0.10' `
            -InterfaceAlias $adapter.InterfaceAlias | Out-String | Out-File -LiteralPath $logFile -Encoding utf8
    }
} catch {
    $_ | Out-String | Out-File -LiteralPath $logFile -Encoding utf8
    exit 1
}
