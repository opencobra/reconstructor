"""Check that Constructor can use the VMH public API (and warm its caches).

    python manage.py check_vmh_api            # smoke lookups + derived InChIKey index
    python manage.py check_vmh_api --no-index # smoke lookups only

Exits with status 1 if any check fails, so it can gate a deployment.
"""
import time

from django.core.management.base import BaseCommand, CommandError

from reactions.utils import vmh_api

WATER_INCHIKEY = 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'


class Command(BaseCommand):
    help = 'Check connectivity to the VMH public API and warm the VMH caches.'

    def add_arguments(self, parser):
        parser.add_argument('--no-index', action='store_true',
                            help='Skip building/loading the derived InChIKey index.')

    def handle(self, *args, **options):
        self.stdout.write(f'VMH API base: {vmh_api.vmh_api_base_url()}')
        self.stdout.write(f'VMH links:    {vmh_api.vmh_public_base_url()}')
        self.stdout.write('VMH API key:  ' + (
            'configured' if vmh_api.vmh_api_key_configured()
            else 'not configured (VMH allows 100 requests/hour per IP without one)'))

        checks = [
            ('metabolite o2', lambda: vmh_api.lookup_metabolite('o2'),
             lambda r: r and r['abbreviation'] == 'o2'),
            ('reaction HEX1', lambda: vmh_api.lookup_reaction('HEX1'),
             lambda r: r and r['abbreviation'] == 'HEX1' and r['formula']),
            ('gene GCK', lambda: vmh_api.lookup_gene('GCK'),
             lambda r: r and r['rows']),
            ('InChIKey of water', lambda: vmh_api.lookup_metabolite_by_inchikey(WATER_INCHIKEY),
             lambda r: r and r['abbreviation']),
        ]
        if options['no_index']:
            pass
        elif not vmh_api.vmh_api_key_configured():
            self.stdout.write('SKIP derived InChIKey index: needs VMH_API_KEY (VMH truncates '
                              'anonymous bulk queries); VMH metabolites without an InChIKey '
                              'are then matched by exact SMILES only')
        else:
            # An empty index is valid (VMH may store InChIKeys for everything
            # one day) but unusual today, so it is reported as a warning below.
            checks.append(('derived InChIKey index', vmh_api.warm_derived_inchikey_index,
                           lambda n: n >= 0))

        failures = 0
        for label, run, ok in checks:
            started = time.time()
            try:
                result = run()
            except vmh_api.VMHUnavailableError as exc:
                failures += 1
                self.stdout.write(self.style.ERROR(f'FAIL {label}: {exc.user_message()}'))
                continue
            elapsed = time.time() - started
            if isinstance(result, int) and result == 0:
                self.stdout.write(self.style.WARNING(
                    f'WARN {label}: 0 keys - VMH exported no metabolites without an '
                    'InChIKey; structure matches for those rely on exact SMILES only'))
            elif ok(result):
                detail = f'{result} keys' if isinstance(result, int) else 'ok'
                self.stdout.write(self.style.SUCCESS(f'OK   {label} ({detail}, {elapsed:.1f}s)'))
            else:
                failures += 1
                self.stdout.write(self.style.ERROR(f'FAIL {label}: unexpected answer {result!r}'))
        if failures:
            raise CommandError(f'{failures} VMH check(s) failed')
