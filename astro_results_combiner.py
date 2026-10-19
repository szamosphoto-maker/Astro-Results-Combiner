"""Astro Results Combiner – iTelescope/Seestar results to VOTable.
Requires: pandas, astropy; tkinter (bundled with regular Windows Python).
"""
from __future__ import annotations

import re
import threading
import queue
import sys
from datetime import datetime, date
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import pandas as pd
from astropy.table import Table

KINDS = {
    'binary': ('binary_results.csv', 'binary_results_osszesitett.vo'),
    'stellar': ('stellar_results.csv', 'stellar_results_osszesitett.vo'),
    'measurement': ('*_measurement.txt', 'measurement_results_osszesitett.vo'),
    'finder': ('*_WDS_Finder.csv', 'WDS_Finder_osszesitett.vo'),
    'rejection': ('*_rejection_diagnostic.txt', 'rejection_diagnostics_osszesitett.vo'),
}
META_COLUMNS = {'Observation_Date', 'Observation_Folder', 'Source_Subfolder',
                'Source_Folder', 'Source_File', 'UTC', 'Object', 'Quality', 'Analyzer', 'FITS_File', 'Rejection', 'Rejection_Type', 'Raw_Diagnostic', 'FITS'}
FORCE_STRING = {'wds', 'wds_id', 'source_id', 'gaia_pri', 'gaia_sec', 'gaia_source_id'}
OPTIONAL_BINARY = ['pmra_1', 'pmdec_1', 'pmra_2', 'pmdec_2']
MEASUREMENT_FIELDS = {
    'Object': 'Object', 'WDS': 'WDS', 'Date': 'UTC', 'PA': 'PA_deg',
    'Separation': 'Separation_arcsec', 'Mag1 (WDS)': 'Mag1_WDS',
    'Mag2 (calculated)': 'Mag2_calculated', 'Delta mag': 'Delta_mag',
    'Quality': 'Quality'
}


def file_date_folder(path: Path, root: Path):
    """Find nearest parent starting with YYYY-MM-DD (including root)."""
    p = path.parent
    while True:
        name = p.name
        if re.match(r'^\d{4}-\d{2}-\d{2}', name):
            try:
                return p, date.fromisoformat(name[:10])
            except ValueError:
                pass
        if p == root or p.parent == p:
            break
        p = p.parent
    return None, None


def clean_column_name(name):
    """Clean BOM artifacts from a CSV column header, preserving real names."""
    name = str(name).strip().lstrip('\ufeff')
    # UTF-8 BOM misread using cp1250 / latin-1, and replacement-char variants.
    for prefix in ('ï»¿', 'ď»ż', 'ď»?', '???', '?»?', '�»�', '���'):
        if name.startswith(prefix):
            name = name[len(prefix):]
    return name.lstrip('\ufeff').strip()



def normalize_identifier_columns(df):
    """Merge duplicated/misdecoded identifier headers without losing values.

    Some older CSV exports contain names like '?»?wds_id' alongside
    the correct 'wds_id'. Drop an alias only after copying any populated
    values into its canonical counterpart; conflicting values are an error.
    """
    df = df.copy()
    for canonical in ('wds_id', 'source_id'):
        aliases = [col for col in df.columns
                   if col != canonical and clean_column_name(col).lower() == canonical]
        # Also recognize anomalous prefixes preceding these exact field names.
        aliases += [col for col in df.columns if col != canonical and col not in aliases
                    and re.fullmatch(r'[^A-Za-z0-9_]+(?:wds_id|source_id)', str(col), re.I)
                    and str(col).lower().endswith(canonical)]
        if not aliases:
            continue
        if canonical not in df.columns:
            df = df.rename(columns={aliases.pop(0): canonical})
        main = df[canonical].fillna('').astype(str)
        for alias in aliases:
            extra = df[alias].fillna('').astype(str)
            both = main.str.strip().ne('') & extra.str.strip().ne('')
            if (both & main.ne(extra)).any():
                raise ValueError(f'Eltérő értékek a(z) {canonical} és {alias} oszlopokban')
            main = main.mask(main.str.strip().eq(''), extra)
            df = df.drop(columns=[alias])
        df[canonical] = main
    return df


def read_csv(path):
    last_error = None
    for encoding in ('utf-8-sig', 'utf-8', 'cp1250', 'latin-1'):
        try:
            # Read as text first: prevent precision loss on 19-digit Gaia IDs.
            df = pd.read_csv(path, encoding=encoding, sep=None, engine='python',
                             dtype=str, keep_default_na=False)
            # Remove BOM and its common encoding-corrupted variants from headers.
            # The values themselves (including 19-digit Gaia IDs) stay untouched.
            df.columns = [clean_column_name(name) for name in df.columns]
            if not df.columns.is_unique:
                raise ValueError('Ismétlődő oszlopnév a fejléc tisztítása után')
            return df, encoding
        except Exception as exc:
            last_error = exc
    raise last_error


def parse_measurement(path):
    last_error = None
    for encoding in ('utf-8-sig', 'utf-8', 'cp1250', 'latin-1'):
        try:
            raw = path.read_text(encoding=encoding)
            break
        except UnicodeError as exc:
            last_error = exc
    else:
        raise last_error or ValueError('Nem olvasható a szövegfájl')

    data = {'Analyzer': '', **{v: '' for v in MEASUREMENT_FIELDS.values()}}
    for line in raw.splitlines():
        line = line.strip()
        if line.startswith('Seestar Double Star Analyzer'):
            data['Analyzer'] = line
        if ':' not in line:
            continue
        key, value = line.split(':', 1)
        key, value = key.strip(), value.strip()
        if key not in MEASUREMENT_FIELDS:
            continue
        field = MEASUREMENT_FIELDS[key]
        if field == 'UTC':
            data[field] = re.sub(r'\s+UTC$', '', value, flags=re.IGNORECASE)
        elif field in {'PA_deg', 'Separation_arcsec', 'Mag1_WDS', 'Mag2_calculated', 'Delta_mag'}:
            match = re.search(r'[-+]?\d+(?:[.,]\d+)?', value)
            data[field] = match.group().replace(',', '.') if match else ''
        else:
            data[field] = re.sub(r'\s+', ' ', value)
    if not all(data.get(k) for k in ('Object', 'WDS', 'UTC', 'PA_deg', 'Separation_arcsec')):
        raise ValueError('Hiányzó kötelező mező: Object, WDS, Date, PA vagy Separation')
    try:
        dt = datetime.fromisoformat(data['UTC'])
    except ValueError as exc:
        raise ValueError('Érvénytelen UTC időpont') from exc
    return data, dt.date(), encoding


def parse_rejection(path: Path):
    """One row per rejected pair; retain original text and detailed diagnostic fields."""
    last_error = None
    for encoding in ('utf-8-sig', 'utf-8', 'cp1250', 'latin-1'):
        try:
            raw = path.read_text(encoding=encoding)
            break
        except UnicodeError as exc:
            last_error = exc
    else:
        raise last_error or ValueError('Nem olvasható a diagnosztikai fájl')

    data = {'Object': '', 'WDS': '', 'FITS_File': '', 'FWHM_px': '',
            'Rejection': '', 'Rejection_Type': '', 'Model_Improvement_pct': '',
            'B_peak_noise': '', 'Raw_Diagnostic': raw}
    for label in ('first', 'second'):
        data[f'{label}_WDS_x_px'] = ''
        data[f'{label}_WDS_y_px'] = ''
        for radius in (3, 6):
            for field in ('peak_x_px', 'peak_y_px', 'method', 'SNR'):
                data[f'{label}_r{radius}_{field}'] = ''

    num = r'([-+]?\d+(?:[.,]\d+)?)'
    active = None
    for line in raw.splitlines():
        stripped = line.strip()
        if stripped.startswith('Object:'):
            data['Object'] = re.sub(r'\s+', ' ', stripped.partition(':')[2].strip())
        elif stripped.startswith('WDS:'):
            data['WDS'] = stripped.partition(':')[2].strip()
        elif stripped.startswith('FITS:'):
            data['FITS_File'] = stripped.partition(':')[2].strip()
        elif stripped.startswith('FWHM estimate:'):
            m = re.search(num, stripped.partition(':')[2])
            if m: data['FWHM_px'] = m.group(1).replace(',', '.')
        elif stripped.startswith('Rejection:'):
            data['Rejection'] = stripped.partition(':')[2].strip()
        elif stripped.startswith(('first:', 'second:')):
            active = stripped.split(':', 1)[0]
            m = re.search(r'WDS predicted pixel=\(\s*'+num+r'\s*,\s*'+num+r'\s*\)', stripped)
            if m:
                data[f'{active}_WDS_x_px'] = m.group(1).replace(',', '.')
                data[f'{active}_WDS_y_px'] = m.group(2).replace(',', '.')
        elif active and stripped.startswith('radius='):
            m = re.search(r'radius=([0-9.]+)px\s+peak=\(\s*'+num+r'\s*,\s*'+num+r'\s*\)\s+method=(.*?)\s+SNR='+num, stripped)
            if m:
                rad = float(m.group(1))
                if rad in (3, 6):
                    tag = f'{active}_r{int(rad)}_'
                    data[tag+'peak_x_px'] = m.group(2).replace(',', '.')
                    data[tag+'peak_y_px'] = m.group(3).replace(',', '.')
                    data[tag+'method'] = m.group(4).strip()
                    data[tag+'SNR'] = m.group(5).replace(',', '.')
    if not data['Object'] or not data['WDS'] or not data['Rejection']:
        raise ValueError('Hiányzó kötelező mező: Object, WDS vagy Rejection')
    data['Rejection_Type'] = data['Rejection'].split(':', 1)[0].strip()
    m = re.search(r'improvement\s+'+num+r'\s*%', data['Rejection'], re.I)
    if m: data['Model_Improvement_pct'] = m.group(1).replace(',', '.')
    m = re.search(r'B peak/noise\s+'+num, data['Rejection'], re.I)
    if m: data['B_peak_noise'] = m.group(1).replace(',', '.')

    # FITS filename timestamp is the available date surrogate, not a measured UTC.
    m = re.search(r'(20\d{6})-\d{6}', data['FITS_File'])
    obs_date = datetime.strptime(m.group(1), '%Y%m%d').date() if m else None
    return data, obs_date, encoding


def finder_details(csv_path: Path):
    """Read the companion Finder report; use its FITS observation date if available."""
    txt_path = csv_path.with_suffix('.txt')
    if not txt_path.is_file():
        return None, None, None
    last_error = None
    for encoding in ('utf-8-sig', 'utf-8', 'cp1250', 'latin-1'):
        try:
            content = txt_path.read_text(encoding=encoding)
            break
        except UnicodeError as exc:
            last_error = exc
    else:
        raise last_error or ValueError('Nem olvasható a Finder TXT')
    match = re.search(r'^Observation \(FITS\):\s*(\d{4}-\d{2}-\d{2})', content, re.MULTILINE)
    observed = date.fromisoformat(match.group(1)) if match else None
    return txt_path, content, observed


def finder_fits_name(csv_path: Path, companion_content: str | None):
    if companion_content:
        match = re.search(r'^FITS:\s*(.+)$', companion_content, re.MULTILINE)
        if match:
            return match.group(1).strip()
    return csv_path.name[:-len('_WDS_Finder.csv')] + '.fit'


def make_votable(df, path):
    df = df.copy()
    for name in df.columns:
        lower = name.lower()
        if (name in META_COLUMNS or lower in FORCE_STRING or 'source_id' in lower
                or lower.startswith('gaia') or lower.startswith('wds')):
            df[name] = df[name].fillna('').astype(str)
            continue
        vals = df[name].fillna('').astype(str)
        nonempty = vals.str.strip().ne('')
        if nonempty.any():
            converted = pd.to_numeric(vals.where(nonempty), errors='coerce')
            if converted[nonempty].notna().all():
                df[name] = converted.astype(float)
                continue
        df[name] = vals
    Table.from_pandas(df, index=False).write(path, format='votable', overwrite=True)


def collect(root: Path, kinds: list[str], recursive: bool, use_dates: bool,
            start: date | None, end: date | None, log):
    outdir = root  # Az eredmények közvetlenül a főkönyvtárba kerülnek.
    totals = []
    for kind in kinds:
        pattern, outname = KINDS[kind]
        paths = sorted(root.rglob(pattern) if recursive else root.glob(pattern))
        # A keresési minták nem illeszkednek a generált összesítő fájlokra.
        audits = []
        parts = []
        for path in paths:
            rel = str(path.relative_to(root))
            row = dict(Type=kind, Status='', Source_File=rel, Observation_Date='',
                       Rows_Read=0, Encoding='', Error='')
            try:
                obs_folder, obs_date = file_date_folder(path, root)
                if kind == 'finder':
                    df, enc = read_csv(path)
                    txt_path, txt_content, report_date = finder_details(path)
                    effective_date = report_date or obs_date
                    df.insert(0, 'FITS_File', finder_fits_name(path, txt_content))
                elif kind == 'rejection':
                    diagnostic, file_obs_date, enc = parse_rejection(path)
                    effective_date = file_obs_date or obs_date
                    df = pd.DataFrame([diagnostic])
                elif kind == 'measurement':
                    measurement, inner_date, enc = parse_measurement(path)
                    effective_date = inner_date
                    df = pd.DataFrame([measurement])
                else:
                    df, enc = read_csv(path)
                    effective_date = obs_date
                    if kind == 'binary':
                        for name in OPTIONAL_BINARY:
                            if name not in df.columns:
                                df[name] = ''
                    # Preserve old iTelescope behavior: a missing dated folder
                    # cannot be included in a date-restricted search.
                row['Encoding'] = enc
                row['Rows_Read'] = len(df)
                if not len(df) and kind != 'finder':
                    row['Status'] = 'EMPTY'
                    audits.append(row)
                    continue
                if use_dates and effective_date is None:
                    row['Status'] = 'NO_DATE'
                    row['Error'] = 'Nem határozható meg a mérés dátuma'
                    audits.append(row)
                    continue
                if use_dates and not (start <= effective_date <= end):
                    row['Status'] = 'OUT_OF_RANGE'
                    row['Observation_Date'] = effective_date.isoformat()
                    audits.append(row)
                    continue
                folder = obs_folder or path.parent
                prefix = {
                    'Observation_Date': effective_date.isoformat() if effective_date else '',
                    'Observation_Folder': folder.name,
                    'Source_Subfolder': str(path.parent.relative_to(folder)),
                    'Source_Folder': str(path.parent.relative_to(root)),
                    'Source_File': rel,
                }
                # Guarantee source columns precede measurement columns.
                for index, (name, value) in enumerate(prefix.items()):
                    if name in df.columns:
                        df[name] = value
                    else:
                        df.insert(index, name, value)
                if len(df):
                    parts.append(df)
                if kind == 'finder' and txt_content is None:
                    log(f'FIGYELEM [finder]: hiányzik a TXT párja: {rel}')
                row['Status'] = 'ACCEPTED'
                row['Observation_Date'] = prefix['Observation_Date']
                log(f'BEOLVASVA [{kind}]: {rel} ({len(df)} sor)')
            except Exception as exc:
                row['Status'] = 'ERROR'
                row['Error'] = str(exc)
                log(f'HIBA [{kind}]: {rel}: {exc}')
            audits.append(row)
        # Nincs külön audit CSV: a feldolgozási napló a GUI-ban látható.
        if audits:
            statuses = {}
            for item in audits:
                status = item['Status']
                statuses[status] = statuses.get(status, 0) + 1
            log(f'Állapotok [{kind}]: ' + ', '.join(f'{k}: {v}' for k, v in sorted(statuses.items())))
        if parts:
            merged = pd.concat(parts, ignore_index=True, sort=False)
            merged = normalize_identifier_columns(merged)
            merged = merged.sort_values(['Observation_Date', 'Source_File'], kind='stable')
            # Az iTelescope-típusú binary/stellar táblák neve csak bekapcsolt
            # dátumszűrésnél kapja meg a tényleges dátumhatárokat.
            if use_dates and kind in ('binary', 'stellar'):
                output_name = f'{Path(outname).stem}_{start.isoformat()}_{end.isoformat()}.vo'
            else:
                output_name = outname
            result_file = outdir / output_name
            make_votable(merged, result_file)
            totals.append(f'{kind}: {len(parts)} fájl, {len(merged)} sor → {output_name}')
        else:
            totals.append(f'{kind}: nincs feldolgozható fájl')
        log(totals[-1])
    return outdir, totals


class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Astro Results Combiner')
        # EXE-ben is megjelenik az ikon az ablak címsorában.
        icon_dir = Path(getattr(sys, '_MEIPASS', Path(__file__).resolve().parent))
        try:
            self.iconbitmap(str(icon_dir / 'Astro_Results_Combiner.ico'))
        except (tk.TclError, OSError):
            pass
        self.geometry('760x710')
        self.minsize(640, 560)
        self.q = queue.Queue()
        self.folder = tk.StringVar()
        self.kinds = {k: tk.BooleanVar(value=k in ('binary', 'measurement')) for k in KINDS}
        self.recursive = tk.BooleanVar(value=True)
        self.use_dates = tk.BooleanVar(value=False)
        year = datetime.now().year
        self.fromdate = tk.StringVar(value=f'{year}-01-01')
        self.todate = tk.StringVar(value=f'{year}-12-31')
        self.build()
        self.after(120, self.poll)

    def build(self):
        main = ttk.Frame(self, padding=16)
        main.pack(fill='both', expand=True)
        ttk.Label(main, text='Astro Results Combiner', font=('Segoe UI', 16, 'bold')).pack(anchor='w')
        ttk.Label(main, text='iTelescope és Seestar mérések egyesítése VOTable formátumba').pack(anchor='w', pady=(0, 15))
        folder_frame = ttk.LabelFrame(main, text='Főkönyvtár', padding=10)
        folder_frame.pack(fill='x', pady=5)
        ttk.Entry(folder_frame, textvariable=self.folder).pack(side='left', fill='x', expand=True)
        ttk.Button(folder_frame, text='Tallózás...', command=self.browse).pack(side='left', padx=(8, 0))
        frame = ttk.LabelFrame(main, text='Beolvasandó fájlok (külön VO-táblákba)', padding=10)
        frame.pack(fill='x', pady=5)
        for kind, (pattern, _) in KINDS.items():
            label = ('DuoStar Finder (*_WDS_Finder.csv)' if kind == 'finder' else
                     'Elutasított mérések (*_rejection_diagnostic.txt)' if kind == 'rejection' else pattern)
            ttk.Checkbutton(frame, text=label, variable=self.kinds[kind]).pack(anchor='w')
        frame = ttk.LabelFrame(main, text='Keresési feltételek', padding=10)
        frame.pack(fill='x', pady=5)
        ttk.Checkbutton(frame, text='Almappák átvizsgálása', variable=self.recursive).grid(row=0, column=0, sticky='w', columnspan=3)
        ttk.Checkbutton(frame, text='Dátumszűrés', variable=self.use_dates).grid(row=1, column=0, sticky='w', columnspan=3)
        ttk.Label(frame, text='Kezdő dátum (ÉÉÉÉ-HH-NN):').grid(row=2, column=0, sticky='w', pady=4)
        ttk.Entry(frame, textvariable=self.fromdate, width=14).grid(row=2, column=1, sticky='w', padx=8)
        ttk.Label(frame, text='Záró dátum:').grid(row=3, column=0, sticky='w', pady=4)
        ttk.Entry(frame, textvariable=self.todate, width=14).grid(row=3, column=1, sticky='w', padx=8)
        self.run_button = ttk.Button(main, text='Összesítés indítása', command=self.start)
        self.run_button.pack(anchor='w', pady=12)
        ttk.Label(main, text='Feldolgozási napló:').pack(anchor='w')
        log_frame = ttk.Frame(main)
        log_frame.pack(fill='both', expand=True)
        self.log_text = tk.Text(log_frame, wrap='word', height=10, state='disabled')
        scroll = ttk.Scrollbar(log_frame, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scroll.set)
        self.log_text.pack(side='left', fill='both', expand=True)
        scroll.pack(side='right', fill='y')

    def browse(self):
        p = filedialog.askdirectory(title='Főkönyvtár kiválasztása')
        if p:
            self.folder.set(p)

    def log(self, message):
        self.q.put(('log', message))

    def start(self):
        root = Path(self.folder.get().strip())
        if not root.is_dir():
            messagebox.showerror('Hiba', 'Válassz egy létező főkönyvtárat.')
            return
        kinds = [name for name, flag in self.kinds.items() if flag.get()]
        if not kinds:
            messagebox.showerror('Hiba', 'Válassz legalább egy fájltípust.')
            return
        dates = self.use_dates.get()
        try:
            begin = date.fromisoformat(self.fromdate.get()) if dates else None
            end = date.fromisoformat(self.todate.get()) if dates else None
            if dates and begin > end:
                raise ValueError('A kezdődátum későbbi a záródátumnál.')
        except ValueError as exc:
            messagebox.showerror('Hibás dátum', str(exc))
            return
        self.run_button.configure(state='disabled')
        recursive = self.recursive.get()
        self.log(f'Főkönyvtár: {root}\nKeresés: {", ".join(kinds)}')
        def worker():
            try:
                output, totals = collect(root, kinds, recursive,
                                         dates, begin, end, self.log)
                self.q.put(('done', f'Elkészült. Kimeneti könyvtár:\n{output}'))
            except Exception as exc:
                self.q.put(('error', str(exc)))
        threading.Thread(target=worker, daemon=True).start()

    def poll(self):
        try:
            while True:
                kind, msg = self.q.get_nowait()
                if kind == 'log':
                    self.log_text.configure(state='normal')
                    self.log_text.insert('end', msg + '\n')
                    self.log_text.see('end')
                    self.log_text.configure(state='disabled')
                elif kind == 'done':
                    self.run_button.configure(state='normal')
                    messagebox.showinfo('Összesítés', msg)
                else:
                    self.run_button.configure(state='normal')
                    messagebox.showerror('Feldolgozási hiba', msg)
        except queue.Empty:
            pass
        self.after(120, self.poll)


if __name__ == '__main__':
    Application().mainloop()
