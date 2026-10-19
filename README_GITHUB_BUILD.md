# Astro Results Combiner — Windows EXE

A program a különféle iTelescope/Seestar mérési adatokból TOPCAT-kompatibilis VOTable-okat készít. A kimenet a kiválasztott főkönyvtárba kerül.

## GitHub feltöltés

A repository főkönyvtárában legyen három fájl:

- `astro_results_combiner.py`
- `Astro_Results_Combiner.ico`
- `requirements.txt`

A workflow fájl pedig **pontosan itt** legyen:

- `.github/workflows/build-windows.yml`

A GitHub webes felületén az **Add file → Create new file** módszerrel az új fájl neve megadható így: `.github/workflows/build-windows.yml`; az útvonal mappáit a GitHub létrehozza.

## Fordítás

A GitHub repositoryban nyisd meg az **Actions** fület, bal oldalt válaszd a **Build Astro Results Combiner for Windows** munkafolyamatot, majd **Run workflow**.

Ha a folyamat zöld pipával lefutott, a futás alján az **Artifacts** résznél töltsd le az `Astro_Results_Combiner_Windows` csomagot. Ebben található a `Astro_Results_Combiner.exe`.

Alternatívaként a fő ág megfelelő fájljainak feltöltésekor a build automatikusan is elindul.

A build a GitHub Windows-gépén fut. Az EXE nem igényel külön Python-telepítést. Első indításnál a Windows SmartScreen figyelmeztetést adhat, mert a fájl nincs kódaláírva.

## Forrásból indítás

`pip install -r requirements.txt`

`python astro_results_combiner.py`
