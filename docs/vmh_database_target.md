# VMH Constructor Database Target

The Constructor add-to-VMH path writes through MATLAB, not Django.

Flow:

1. Django prepares JSON files for the selected reactions.
2. Django calls the MATLAB HTTP service.
3. MATLAB runs `updateVMHFromConstructor`.
4. `updateVMHFromConstructor` calls `initialiseMySqlCommand`.
5. `/matlab/toolboxes/vmh_revamped/vmh_db_update/src/initialiseMySqlCommand.m` builds the MySQL command from environment variables.

The target database is controlled by:

```bash
VMH_MYSQL_HOST=127.0.0.1
VMH_MYSQL_PORT=3306
VMH_MYSQL_DATABASE=reconDBtest
VMH_MYSQL_USER=saleh
VMH_MYSQL_PASSWORD=
```

These values live in `.env` for the deployed app. Docker only passes them through.

Recommended database setup:

1. Keep the original VMH database unchanged, for example `reconDB_original`.
2. Create a writable clone for Constructor additions, for example `reconDB_constructor`.
3. Set `VMH_MYSQL_DATABASE=reconDB_constructor` in `.env`.
4. Recreate the `web` and `matlab` services so both see the same target.

Important: the current UI availability checks use the public VMH API. If the constructor clone diverges from public VMH, those checks will not see prior constructor-only additions unless we add a target-database availability check.
