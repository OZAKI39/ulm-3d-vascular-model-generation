"""Append-only, evidence-bound extensions to the original shared GPU pool.

No code path creates approval. Records are accepted only through the separate
explicit registration command after a real user message authorizes this scope.
The original limit and every old charge remain intact; extensions are spent
first, only by their named campaign/task, and never grant automatic retries.
"""
import math
from pathlib import Path
from .common import read_json, sha256_file


def authorization_records(pool):
    records = []
    seen = set()
    for reference in pool.get('authorization_records', []):
        path = Path(reference['path'])
        if sha256_file(path) != reference['sha256']:
            raise ValueError('AUTHORIZATION_RECORD_HASH_MISMATCH')
        record = read_json(path)
        evidence = Path(record['user_message_file'])
        if (record['id'] in seen or record.get('authority') != 'explicit_user_message'
                or sha256_file(evidence) != record['user_message_sha256']
                or evidence.read_text().strip() != record['user_message'].strip()
                or not record['user_message'].strip()):
            raise ValueError('INVALID_USER_AUTHORIZATION_EVIDENCE')
        seconds = record['additional_seconds']
        if isinstance(seconds, bool) or not math.isfinite(seconds) or seconds <= 0:
            raise ValueError('INVALID_ADDITIONAL_SECONDS')
        if record['single_task_limit_s'] != 600 or record['gpu_concurrency'] != 1:
            raise ValueError('UNSUPPORTED_AUTHORIZATION_SCOPE')
        request = read_json(Path(record['budget_request_file']))
        if sha256_file(Path(record['budget_request_file'])) != record['budget_request_sha256']:
            raise ValueError('AUTHORIZATION_REQUEST_HASH_MISMATCH')
        if record['scope'] != request['authorization_scope']:
            raise ValueError('AUTHORIZATION_SCOPE_MISMATCH')
        seen.add(record['id'])
        records.append(record)
    return records


def account_extensions(pool, attempts, campaign, task_id):
    """All charges/reservations remain charged, including failed/crashed jobs."""
    records = authorization_records(pool)
    scopes = {}
    for record in records:
        scope = record['scope']
        key = (str(Path(scope['campaign_directory']).resolve()), scope['task_id'])
        scopes[key] = scopes.get(key, 0.) + record['additional_seconds']
    total = sum(charge for _, _, charge in attempts)
    extra_used = 0.
    eligible_remaining = 0.
    rows = []
    for (directory, task), authorized in scopes.items():
        spent = sum(charge for path, name, charge in attempts if path == directory and name == task)
        debit = min(spent, authorized)
        extra_used += debit
        remaining = max(0., authorized - debit)
        if (directory, task) == (str(Path(campaign).resolve()), task_id):
            eligible_remaining += remaining
        rows.append({'campaign_directory': directory, 'task_id': task,
                     'authorized_s': authorized, 'charged_or_reserved_to_extension_s': debit,
                     'remaining_s': remaining})
    base_used = total - extra_used
    base_remaining = max(0., pool['limit_s'] - base_used)
    global_limit = pool['limit_s'] + sum(scopes.values())
    return {'original_limit_s': pool['limit_s'], 'extra_authorized_gpu_seconds': sum(scopes.values()),
            'original_charged_or_reserved_s': base_used, 'original_remaining_s': base_remaining,
            'extension_accounts': rows, 'limit_s': global_limit,
            'global_remaining_s': max(0., global_limit - total),
            'remaining_s': min(max(0., global_limit - total), base_remaining + eligible_remaining),
            'authorization_records': records}
