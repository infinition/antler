# Intégrité de l’archive de travail

`FILE_MANIFEST.txt`, `MANIFEST.sha256` et `SHA256SUMS.txt` sont régénérés à
partir de tous les fichiers ordinaires de l’archive, sauf :

- les trois fichiers de manifeste eux-mêmes ;
- les caches Python `__pycache__` et `*.pyc` ;
- les logs bruts potentiellement réécrits sous `results/raw/phase4_7/`.

Les JSON sérialisés, y compris les campagnes Phase 4.7 terminées, sont inclus.
Les logs bruts restent hors manifeste tant qu'ils peuvent être prolongés.

Pour vérifier la copie courante :

```powershell
Get-Content -Encoding UTF8 MANIFEST.sha256 | ForEach-Object {
  $hash, $path = $_ -split '\s+', 2
  if ((Get-FileHash -Algorithm SHA256 $path.Substring(2)).Hash.ToLower() -ne $hash) {
    throw "checksum mismatch: $path"
  }
}
```
