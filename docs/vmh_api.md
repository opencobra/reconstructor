# VMH integration

Constructor reads metabolites, reactions and genes from the public VMH API at
`https://www.vmh.life/api/public` (documentation: <https://www.vmh.life/apiDocs>).
VMH retired its legacy `/_api/...` endpoints and the interim `vmh.chatimd.org`
host; see <https://www.vmh.life/apiMigration>. All access goes through
`curationTool/reactions/utils/vmh_api.py`.

## Configuration (`.env`)

| Variable | Default | Purpose |
|---|---|---|
| `VMH_API_KEY` | *(empty)* | Optional VMH `public_unlimited` key, sent as `x-api-key`. Without it VMH allows **100 requests per hour per server IP**, shared by all endpoints. |
| `VMH_BASE_URL` | `https://www.vmh.life` | Site used for links (`/metabolite/<abbr>`, `/reaction/<abbr>`). |
| `VMH_API_BASE_URL` | `https://www.vmh.life/api/public` | API base URL. |
| `VMH_API_TIMEOUT` | `20` | Seconds per request (structure fallbacks use longer timeouts). |
| `VMH_CACHE_TTL` / `VMH_NEGATIVE_CACHE_TTL` | `86400` / `3600` | Seconds to cache answers / "not found" answers. |
| `DJANGO_CACHE_DIR` | `<tmp>/constructor-cache` | Shared file cache for all gunicorn workers (VMH answers, API throttling). |
| `VMH_DERIVED_INDEX` | `true` | Derive InChIKeys for VMH metabolites stored without one (needs `VMH_API_KEY`). |

A value pointing at the retired `vmh.chatimd.org` host is ignored (with a warning).

## Behaviour

- **"Not found" and "VMH unreachable" are different answers.** When VMH cannot
  be asked (timeout, HTTP 429/5xx, unexpected response) curators see the
  reason, e.g. *Could not verify `o2` in VMH: VMH's request limit has been
  reached. Please try again in about 12 minutes.*; the API answers
  `502 upstream_unavailable`. Checks that gate submissions to VMH (free names
  and abbreviations) never treat an unanswered question as "free": nothing is
  added and the reason is shown.
- **Rate limit.** After a 429 every VMH call fails fast until VMH's advertised
  reset time. Answers are cached, so repeated lookups cost no request.
- **Structure matching** (ChEBI, PubChem, SMILES, drawings, MOL files) uses
  VMH's InChIKeys via its SPARQL endpoint, then exact SMILES. About half of
  VMH's metabolites have a SMILES but no InChIKey; with `VMH_API_KEY` set,
  Constructor derives their InChIKeys with RDKit into an index that is built on
  first use (~25 s), saved in `DJANGO_CACHE_DIR/vmh` and refreshed daily in the
  background. Without a key this index is skipped (VMH truncates anonymous bulk
  queries), so those metabolites are matched by exact SMILES only.
- **Abbreviations.** Former (renamed) metabolite abbreviations resolve to the
  current one in verification; VMH matches reaction abbreviations
  case-sensitively, so the upper- and lower-case spellings are tried too.

ChEBI and Rhea are used through their current REST APIs
(`reactions/utils/chebi_api.py`, `reactions/utils/get_from_rhea.py`); the ChEBI
SOAP web service was retired by EBI. ChEBI answers each request in ~15 s, so
compounds are fetched in batches and cached.

## Checking a deployment

```bash
python manage.py check_vmh_api          # inside the web container: docker compose exec web ...
```

It prints the configuration, looks up `o2`, `HEX1`, `GCK` and water's InChIKey,
builds (or loads) the derived index, and exits non-zero if anything fails.
