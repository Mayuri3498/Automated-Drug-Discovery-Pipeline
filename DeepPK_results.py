import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

def ensure_playwright():
    try:
        from playwright.sync_api import sync_playwright
        return sync_playwright
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "playwright"])
        from playwright.sync_api import sync_playwright
        return sync_playwright

def find_chrome():
    candidates = [os.environ.get("CHROME_PATH"), r"C:\Program Files\Google\Chrome\Application\chrome.exe",r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
        os.path.expandvars(r"%LOCALAPPDATA%\Google\Chrome\Application\chrome.exe")]
    for p in candidates:
        if p and Path(p).exists():
            return p
    return None

def norm(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())

def read_top10(path):
    import pandas as pd
    df = pd.read_csv(path)
    if df.empty:
        raise RuntimeError("top10_screened.csv is empty.")

    preferred = ["Molecule", "Molecule Name", "Compound", "Compound Name", "Name", "ID", "Compound ID", "Molecule ID", "Filename", "File", "SDF", "Ligand"]
    normalized = {norm(c): c for c in df.columns}
    col = None

    for name in preferred:
        if norm(name) in normalized:
            col = normalized[norm(name)]
            break

    if col is None:
        col = df.columns[0]

    values = [str(x).strip() for x in df[col].dropna().tolist()]

    return values

def load_sdfs(folder):
    from rdkit import Chem
    records = []

    for p in sorted(Path(folder).glob("*.sdf")):
        supplier = Chem.SDMolSupplier(str(p), removeHs=False)

        for index, mol in enumerate(supplier):
            if mol is None:
                continue

            props = {}
            for key in mol.GetPropNames():
                try:
                    props[key] = mol.GetProp(key)
                except Exception:
                    pass

            records.append({"path": p, "stem": p.stem, "record_index": index, "smiles": Chem.MolToSmiles(mol, canonical=True), "props": props})

    if not records:
        raise RuntimeError(f"No readable SDF molecules found in {folder}.")

    return records

def match_records(top_values, records):
    lookup = {}
    for record in records:
        values = {norm(record["stem"]), norm(record["path"].name), norm(record["path"].stem)}

        for value in record["props"].values():
            if isinstance(value, str):
                values.add(norm(value))

        for value in values:
            if value:
                lookup.setdefault(value, []).append(record)

    selected = []
    used = set()

    for value in top_values:
        candidates = lookup.get(norm(value), [])
        if not candidates:
            stripped = re.sub(r"\.(sdf|sd)$", "", value, flags=re.I)
            candidates = lookup.get(norm(stripped), [])

        if not candidates:
            raise RuntimeError(f"Could not match '{value}' to an SDF molecule in the molecules folder.")

        record = next((item for item in candidates if id(item) not in used), candidates[0])
        selected.append(record)
        used.add(id(record))

    return selected

def write_smiles(selected, output):
    path = output / "molecules_smiles.smi"
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("SMILES\n")
        for record in selected:
            handle.write(record["smiles"] + "\n")

    return path

def build_metadata(selected):
    rows = []
    for order, record in enumerate(selected, 1):
        row = {"Input_Order": order, "SDF_File": record["path"].name, "SDF_Record": record["record_index"], "Canonical_SMILES": record["smiles"]}
        for key, value in record["props"].items():
            if key not in row:
                row[f"SDF_{key}"] = value

        rows.append(row)
    return rows

def find_result_url(text):
    pattern = r"https?://biosig\.lab\.uq\.edu\.au/deeppk/prediction_results/admet_[A-Za-z0-9_.-]+"
    match = re.search(pattern, text)
    return match.group(0) if match else None

def page_has_result(page):
    try:
        if page.locator("table").count() > 0:
            for i in range(page.locator("table").count()):
                try:
                    if page.locator("table").nth(i).locator("tr").count() >= 2:
                        return True
                except Exception:
                    pass
    except Exception:
        pass

    try:
        for locator in ["a[href*='.csv']", "a[href*='download']", "a:has-text('CSV')", "a:has-text('Download')", "button:has-text('CSV')", "button:has-text('Download')"]:
            if page.locator(locator).count() > 0:
                return True
    except Exception:
        pass

    return False

def safe_revisit(page, url):
    try:
        if page.is_closed():
            return False

        page.goto(url, wait_until="commit")

        try:
            page.wait_for_load_state("domcontentloaded")
        except Exception as exc:
            message = str(exc)
            if "ERR_ABORTED" not in message and "frame was detached" not in message:
                print(f"Result-page load warning: {exc}")
        return "/prediction_results/" in page.url

    except Exception as exc:
        message = str(exc)

        if ("ERR_ABORTED" in message or "frame was detached" in message or "page is navigating" in message):
            print("Result-page navigation transition detected; continuing.")
            return True

        print(f"Result-page revisit warning: {exc}")
        return False

def wait_for_results(page, result_url):
    last_error = None
    revisit_interval = 15
    last_revisit = 0

    while True:
        try:
            if page.is_closed():
                raise RuntimeError("The Deep-PK page was closed unexpectedly.")

            current = page.url

            if "/prediction_results/" not in current:
                safe_revisit(page, result_url)

            try:
                page.wait_for_load_state("domcontentloaded")
            except Exception as exc:
                message = str(exc)
                if "ERR_ABORTED" not in message and "frame was detached" not in message:
                    print(f"DOM wait warning: {exc}")

            if page.is_closed():
                raise RuntimeError("The Deep-PK page was closed unexpectedly.")
            current = page.url

            if "/prediction_results/" not in current:
                page.wait_for_timeout(2000)
                continue

            if page_has_result(page):
                print("Result content detected.")
                return True

            try:
                body = page.locator("body").inner_text().lower()
            except Exception as exc:
                message = str(exc)
                if "navigating" in message or "detached" in message or "ERR_ABORTED" in message:
                    print("Page is transitioning; waiting before inspection...")
                    page.wait_for_timeout(3000)
                    continue
                raise

            if "internal server error" in body:
                raise RuntimeError("Deep-PK returned an Internal Server Error page.")

            if "bad request" in body:
                raise RuntimeError("Deep-PK returned a Bad Request page.")

            if "done" in body:
                print("Deep-PK reports the job is done. Looking for the result content...")

            print("Waiting for Deep-PK results...")
            now = time.time()

            if now - last_revisit >= revisit_interval:
                safe_revisit(page, result_url)
                last_revisit = now
            page.wait_for_timeout(5000)

        except Exception as exc:
            last_error = exc
            message = str(exc)

            if ("ERR_ABORTED" in message or "frame was detached" in message or "page is navigating" in message or "Unable to retrieve content" in message):
                print("Browser navigation transition detected; recovering...")
                page.wait_for_timeout(3000)
                safe_revisit(page, result_url)
                page.wait_for_timeout(3000)
                continue

            raise

def extract_result_csv(page, output, metadata):
    import pandas as pd
    from urllib.parse import urljoin
    output.mkdir(parents=True, exist_ok=True)

    data = page.evaluate("""
() => {
    const result = [];
    if (!window.jQuery) return result;
    const tables = document.querySelectorAll('table');
    for (const table of tables) {
        try {
            if (window.jQuery.fn.dataTable && window.jQuery.fn.dataTable.isDataTable(table)) {
                const dt = window.jQuery(table).DataTable();
                const info = dt.page.info();
                const headers = Array.from(table.querySelectorAll('thead th')).map(x => x.innerText.trim());
                const rows = dt.rows({search: 'none'}).data().toArray();
                result.push({
                    id: table.id || '',
                    headers: headers,
                    recordsTotal: info.recordsTotal,
                    rows: rows
                });
            }
        } catch (e) {}
    }
    return result;
}
""")

    if data:
        candidates = [x for x in data if x.get('rows')]
        if candidates:
            best = max(candidates, key=lambda x: len(x.get('rows', [])))
            rows = best.get('rows', [])
            headers = best.get('headers', [])
            records_total = int(best.get('recordsTotal') or 0)

            if records_total and len(rows) != records_total:
                raise RuntimeError(f"Deep-PK reports {records_total} result rows but DataTables exposed {len(rows)} rows.")
            normalized_rows = []
            max_cols = max(len(r) if isinstance(r, list) else len(r.keys()) for r in rows)
            if len(headers) < max_cols:
                headers = headers + [f"Column_{i+1}" for i in range(len(headers), max_cols)]
            headers = headers[:max_cols]

            for row in rows:
                if isinstance(row, dict):
                    values = list(row.values())
                else:
                    values = list(row)
                values += [''] * (max_cols - len(values))
                normalized_rows.append(values[:max_cols])

            table = pd.DataFrame(normalized_rows, columns=headers)
            if not table.empty:
                raw = output / "DeepPK_results_raw.csv"
                table.to_csv(raw, index=False)
                print(f"Extracted complete DataTables result: {len(table)} rows x {len(table.columns)} columns")
                return raw

    html = page.content()
    links = []
    for a in page.locator("a").all():
        try:
            text = (a.inner_text() or "").strip()
            href = a.get_attribute("href") or ""
            onclick = a.get_attribute("onclick") or ""
            links.append((text, href, onclick))
        except Exception:
            pass

    for text, href, onclick in links:
        combined = f"{text} {href} {onclick}".lower()
        if ".csv" in combined:
            target = urljoin(page.url, href) if href else None
            if target:
                try:
                    response = page.request.get(target)
                    if response.ok and len(response.body()) > 50:
                        raw = output / "DeepPK_results_raw.csv"
                        raw.write_bytes(response.body())
                        df = pd.read_csv(raw)
                        if not df.empty:
                            print(f"Downloaded result CSV: {len(df)} rows x {len(df.columns)} columns")
                            return raw
                except Exception:
                    pass

    debug = output / "DeepPK_results_debug.html"
    debug.write_text(html, encoding="utf-8")
    raise RuntimeError(
        "Deep-PK result page was reached, but its complete DataTables result data could not be extracted. "
        f"Debug page saved to: {debug}"
    )

def extract_all_focus_links(page):
    rows = page.evaluate("""
() => {
    const result = [];
    if (!window.jQuery) return result;
    for (const table of document.querySelectorAll("table")) {
        try {
            if (window.jQuery.fn.dataTable && window.jQuery.fn.dataTable.isDataTable(table)) {
                const dt = window.jQuery(table).DataTable();
                const data = dt.rows().data().toArray();
                result.push({
                    recordsTotal: dt.page.info().recordsTotal,
                    data: data
                });
            }
        } catch (e) {}
    }
    return result;
}
""")
    if not rows:
        raise RuntimeError("Deep-PK result DataTable was not found.")

    best = max(rows, key=lambda x: len(x.get("data", [])))
    records_total = int(best.get("recordsTotal") or 0)
    data = best.get("data") or []

    if records_total != len(data):
        raise RuntimeError(f"Deep-PK reports {records_total} result rows but the DataTable returned {len(data)} rows.")
    links = []
    for row_index, row in enumerate(data, 1):
        cells = row if isinstance(row, list) else list(row.values())
        focus_url = None
        for cell in cells:
            html = str(cell)
            match = re.search(r"""href=["\']([^"\']*horizontal_view[^"\']*)["\']""", html, re.I)
            if match:
                focus_url = match.group(1)
                break
        if focus_url:
            if focus_url.startswith("http://") or focus_url.startswith("https://"):
                url = focus_url
            else:
                url = "https://biosig.lab.uq.edu.au" + (focus_url if focus_url.startswith("/") else "/" + focus_url)
            links.append((row_index, url))

    if len(links) != records_total:
        raise RuntimeError(f"Deep-PK returned {records_total} molecules but only {len(links)} PK Focus links were found.")

    return links

def clean_text(value):
    return re.sub(r"\s+", " ", str(value).strip())

def unique_column_name(name, used):
    base = clean_text(name)
    if not base:
        base = "Unnamed_Property"
    candidate = base
    n = 2
    while candidate in used:
        candidate = f"{base}_{n}"
        n += 1
    used.add(candidate)
    return candidate

def find_section_table(soup, heading_text, header_keywords=None):
    headings = soup.find_all(re.compile(r"^h[1-6]$"))
    target = None
    wanted = clean_text(heading_text).lower()
    for h in headings:
        if clean_text(h.get_text(" ", strip=True)).lower() == wanted:
            target = h
            break
    if target is None:
        return None
    table = target.find_next("table")
    if table is None:
        return None
    if header_keywords:
        headers = [clean_text(x.get_text(" ", strip=True)).lower() for x in table.find_all("th")]
        joined = " | ".join(headers)
        if not all(k.lower() in joined for k in header_keywords):
            candidate = table
            for t in target.find_all_next("table"):
                headers2 = [clean_text(x.get_text(" ", strip=True)).lower() for x in t.find_all("th")]
                joined2 = " | ".join(headers2)
                if all(k.lower() in joined2 for k in header_keywords):
                    candidate = t
                    break
            table = candidate
    return table

def parse_molecule_properties_table(table):
    properties = {}
    if table is None:
        return properties
    rows = table.find_all("tr")
    for tr in rows:
        cells = [clean_text(c.get_text(" ", strip=True)) for c in tr.find_all(["th", "td"])]
        if len(cells) < 2:
            continue
        key = cells[0]
        value = cells[1]
        if key.lower() in ("descriptor", "value"):
            continue
        key = re.sub(r"\s*\([^)]*\)\s*$", "", key).strip()
        if key:
            properties[key] = value
    return properties

def parse_predicted_properties_table(table):
    detailed = []
    if table is None:
        return detailed
    rows = table.find_all("tr")
    headers = []
    for tr in rows:
        ths = tr.find_all("th")
        if ths:
            headers = [clean_text(x.get_text(" ", strip=True)) for x in ths]
            break

    for tr in rows:
        cells = [clean_text(c.get_text(" ", strip=True)) for c in tr.find_all(["th", "td"])]
        if not cells or (headers and cells == headers):
            continue
        if len(cells) < 4:
            continue
        if headers and len(cells) < len(headers):
            cells += [""] * (len(headers) - len(cells))
        if len(cells) >= 7:
            category = cells[0]
            property_name = cells[1]
            predicted_value = cells[2]
            unit = cells[3]
            confidence = cells[4]
            pred_range = cells[5]
            mw_range = cells[6]
        else:
            category = cells[0]
            property_name = cells[1]
            predicted_value = cells[2]
            unit = cells[3] if len(cells) > 3 else ""
            confidence = cells[4] if len(cells) > 4 else ""
            pred_range = cells[5] if len(cells) > 5 else ""
            mw_range = cells[6] if len(cells) > 6 else ""
        if property_name.lower() in ("property name", ""):
            continue
        detailed.append({"PK_and_Toxicity_Category": category, "Property_Name": property_name, "Predicted_Value": predicted_value, "Property_Unit": unit,
            "Predictive_Confidence": confidence, "Prediction_Within_Training_Range": pred_range, "MW_Within_Training_Range": mw_range})
    return detailed

def extract_focus_information(page, debug_dir=None, molecule_number=None):
    from bs4 import BeautifulSoup
    try:
        page.wait_for_function(
            """() => document.body && (
                document.body.innerText.includes('Molecule Properties') ||
                document.body.innerText.includes('Predicted Molecule Properties')
            )""",
        )
    except Exception:
        pass

    while True:
        try:
            body_text = page.locator("body").inner_text()
            if "Predicted Molecule Properties" in body_text and "Molecule Properties" in body_text:
                html = page.content()
                soup = BeautifulSoup(html, "html.parser")

                molecule_table = find_section_table(soup, "Molecule Properties", header_keywords=["Descriptor", "Value"])
                predicted_table = find_section_table(soup, "Predicted Molecule Properties", header_keywords=["PK and Toxicity", "Property Name", "Predicted Value"])
                molecule_props = parse_molecule_properties_table(molecule_table)
                predicted_props = parse_predicted_properties_table(predicted_table)

                if molecule_props and predicted_props:
                    return molecule_props, predicted_props

                all_tables = soup.find_all("table")
                candidate_molecule = None
                candidate_predicted = None

                for table in all_tables:
                    headers = [clean_text(x.get_text(" ", strip=True)).lower() for x in table.find_all("th")]
                    joined = " | ".join(headers)
                    if "descriptor" in joined and "value" in joined and candidate_molecule is None:
                        candidate_molecule = table
                    if "property name" in joined and "predicted value" in joined and candidate_predicted is None:
                        candidate_predicted = table

                if candidate_molecule is not None and candidate_predicted is not None:
                    molecule_props = parse_molecule_properties_table(candidate_molecule)
                    predicted_props = parse_predicted_properties_table(candidate_predicted)
                    if molecule_props and predicted_props:
                        return molecule_props, predicted_props
        except Exception:
            pass

        page.wait_for_timeout(1000)

    html = ""
    try:
        html = page.content()
    except Exception:
        pass

    if debug_dir:
        debug_dir = Path(debug_dir)
        debug_dir.mkdir(parents=True, exist_ok=True)
        suffix = f"_{molecule_number}" if molecule_number is not None else ""
        (debug_dir / f"focus_failed{suffix}.html").write_text(html, encoding="utf-8")
        try:
            page.screenshot(path=str(debug_dir / f"focus_failed{suffix}.png"), full_page=True)
        except Exception:
            pass

    return {}, []

def flatten_focus_result(number, label, smiles, molecule_props, predicted_props):
    row = {"Molecule_Number": number, "Molecule_Label": label, "Canonical_SMILES": smiles, "Result_Type": "PK Focus"}
    descriptor_map = {"Molecular Weight": "MW", "LogP": "LogP", "#Rotatable Bonds": "Rotatable Bonds", "#Acceptors": "Acceptors", "#Donors": "Donors", "Surface Area": "SA"}
    for source_key, output_key in descriptor_map.items():
        if source_key in molecule_props:
            row[output_key] = molecule_props[source_key]
        else:
            alternate = next((k for k in molecule_props if clean_text(k).lower() == source_key.lower()), None)
            if alternate:
                row[output_key] = molecule_props[alternate]
            else:
                row[output_key] = ""

    used = set(row.keys())

    for item in predicted_props:
        property_name = clean_text(item["Property_Name"])
        if not property_name:
            continue
        column = unique_column_name(property_name, used)
        row[column] = item["Predicted_Value"]

    return row

def save_focus_flat(rows, output_file):
    import pandas as pd
    if not rows:
        pd.DataFrame().to_csv(output_file, index=False)
        return 0
    df = pd.DataFrame(rows)
    preferred = ["Molecule_Number", "Molecule_Label", "Canonical_SMILES", "Result_Type", "MW", "LogP", "Rotatable Bonds", "Acceptors", "Donors", "SA"]
    ordered = [c for c in preferred if c in df.columns]
    ordered += [c for c in df.columns if c not in ordered]
    df = df[ordered]
    df.to_csv(output_file, index=False)
    return len(df)

def save_focus_detailed(rows, output_file):
    import pandas as pd
    if not rows:
        pd.DataFrame().to_csv(output_file, index=False)
        return 0
    df = pd.DataFrame(rows)
    df.to_csv(output_file, index=False)
    return len(df)

def collect_focus_for_all_molecules(page, output, selected, metadata):
    focus_links = extract_all_focus_links(page)

    if len(focus_links) != len(selected):
        raise RuntimeError(f"STOP: expected {len(selected)} PK Focus links but detected {len(focus_links)}.")

    with (output / "DeepPK_focus_links.json").open("w", encoding="utf-8") as handle:
        json.dump(
            [
                {
                    "Molecule_Number": i,
                    "SDF_File": metadata[i - 1]["SDF_File"],
                    "Canonical_SMILES": metadata[i - 1]["Canonical_SMILES"],
                    "Focus_URL": url
                }
                for i, url in focus_links
            ],
            handle,
            indent=2,
            ensure_ascii=False
        )

    flat_rows = []
    detailed_rows = []
    failures = []

    for position, focus_url in focus_links:
        meta = metadata[position - 1]
        label = meta["SDF_File"]
        smiles = meta["Canonical_SMILES"]

        print(f"[PK Focus {position}/{len(selected)}] {label}")

        try:
            page.goto(focus_url, wait_until="domcontentloaded",)
            molecule_props, predicted_props = extract_focus_information(page, debug_dir=output / "PK_Focus_debug", molecule_number=position)
            if not molecule_props:
                raise RuntimeError("Molecule Properties section was not detected.")
            if not predicted_props:
                raise RuntimeError("Predicted Molecule Properties section was not detected.")

            flat_rows.append(flatten_focus_result(position, label, smiles, molecule_props, predicted_props))
            for item in predicted_props:
                detailed_rows.append({"Molecule_Number": position, "Molecule_Label": label, "Canonical_SMILES": smiles, **item})

        except Exception as exc:
            failures.append({"Molecule_Number": position, "SDF_File": label, "Canonical_SMILES": smiles, "Focus_URL": focus_url, "Error": str(exc)})
            print(f"  PK Focus failed: {exc}")

        if position % 10 == 0 or position == len(selected):
            save_focus_flat(flat_rows, output / "PK_Focus_all_molecules.csv")
            save_focus_detailed(detailed_rows, output / "PK_Focus_detailed_all_molecules.csv")
            with (output / "DeepPK_focus_failures.json").open("w", encoding="utf-8") as handle:
                json.dump(failures, handle, indent=2, ensure_ascii=False)

    flat_count = save_focus_flat(flat_rows, output / "PK_Focus_all_molecules.csv")
    detailed_count = save_focus_detailed(detailed_rows, output / "PK_Focus_detailed_all_molecules.csv")

    if failures:
        raise RuntimeError(f"PK Focus extraction finished with {len(failures)} molecule failures. See DeepPK_focus_failures.json.")

    if flat_count != len(selected):
        raise RuntimeError(f"PK Focus verification failed: expected {len(selected)} molecule rows, got {flat_count}.")

    return flat_count, detailed_count

def combine_results(raw_csv, metadata, output):
    import pandas as pd
    results = pd.read_csv(raw_csv)
    meta = pd.DataFrame(metadata)
    smiles_col = None

    for column in results.columns:
        if "smiles" in norm(column):
            smiles_col = column
            break

    if smiles_col and not results.empty:
        results["_smiles_key"] = results[smiles_col].astype(str).map(norm)
        meta["_smiles_key"] = meta["Canonical_SMILES"].astype(str).map(norm)
        combined = results.merge(meta, on="_smiles_key", how="left")
        combined.drop(columns=["_smiles_key"], inplace=True)

    else:
        combined = results.copy()

        for column in meta.columns:
            if column not in combined.columns:
                values = meta[column].tolist()
                combined[column] = values[:len(combined)]

    combined.to_csv(output, index=False)
    return output

def main():
    global result_dir_global
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdf_folder", required=True)
    parser.add_argument("--top10_csv", required=True)
    parser.add_argument("--output_folder", required=True)

    args = parser.parse_args()
    output = Path(args.output_folder)
    output.mkdir(parents=True, exist_ok=True)
    result_dir_global = str(output)

    print("DeepPK_AUTOMATED_BROWSER_ADMET_PK_FOCUS")
    print(f"Preparing all molecules from {args.top10_csv}...")

    top_values = read_top10(args.top10_csv)
    records = load_sdfs(args.sdf_folder)
    selected = match_records(top_values, records)
    smiles_file = write_smiles(selected, output)
    metadata = build_metadata(selected)

    with (output / "selected_molecules.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
    print(f"SMILES file created: {smiles_file}")

    sync_playwright = ensure_playwright()
    chrome = find_chrome()

    with sync_playwright() as playwright:
        launch_options = {"headless": False, "args": ["--disable-blink-features=AutomationControlled"]}
        if chrome:
            launch_options["executable_path"] = chrome

        browser = playwright.chromium.launch(**launch_options)

        context = browser.new_context(accept_downloads=True, viewport={"width": 1400, "height": 900})

        page = context.new_page()

        page.set_default_timeout(0)
        page.set_default_navigation_timeout(0)

        print("Opening Deep-PK...")
        page.goto("https://biosig.lab.uq.edu.au/deeppk/prediction", wait_until="domcontentloaded",)

        file_input = page.locator('input[name="smiles_file"]')
        file_input.wait_for(state="attached")
        file_input.set_input_files(str(smiles_file))
        admet_button = page.locator('button[name="pred_type"][value="admet"]')
        admet_button.wait_for(state="visible")
        print("Submitting ADMET prediction...")

        try:
            admet_button.click(no_wait_after=True)
        except Exception:
            admet_button.click()

        start = time.time()
        result_url = None

        try:
            page.wait_for_url("**/prediction_results/**", wait_until="commit")
            result_url = page.url
        except Exception as exc:
            print(f"Result navigation warning: {exc}")

        if result_url is None or "/prediction_results/" not in result_url:
            while result_url is None:
                try:
                    current = page.url

                    if "/prediction_results/" in current:
                        result_url = current
                        break

                    if page.is_closed():
                        raise RuntimeError("The Deep-PK browser page was closed unexpectedly.")

                    try:
                        page.wait_for_load_state("domcontentloaded")
                    except Exception:
                        pass

                    current = page.url

                    if "/prediction_results/" in current:
                        result_url = current
                        break

                    print("Waiting for Deep-PK to create the result page...")
                    page.wait_for_timeout(2000)

                except Exception as exc:
                    message = str(exc)

                    if ("ERR_ABORTED" in message or "frame was detached" in message or "page is navigating" in message or "Target page" in message):
                        print("Browser navigation is still changing; waiting safely...")
                        page.wait_for_timeout(3000)
                        continue

                    raise

        if not result_url:
            debug = output / "DeepPK_submission_debug.html"

            try:
                if not page.is_closed():
                    try:
                        page.wait_for_load_state("domcontentloaded")
                    except Exception:
                        pass
                    debug.write_text(page.content(), encoding="utf-8")
            except Exception:
                debug.write_text("Unable to capture page HTML because navigation was still active.", encoding="utf-8")

            raise RuntimeError(
                "Deep-PK submission did not produce a prediction_results URL. "
                f"Debug information saved to: {debug}"
            )

        print(f"Deep-PK job submitted successfully:\n{result_url}")

        (output / "DeepPK_result_url.txt").write_text(result_url + "\n", encoding="utf-8")

        wait_for_results(page,result_url)

        raw_csv = extract_result_csv(page, output, metadata)
        final_csv = output / "DeepPK_all_properties.csv"

        combine_results(raw_csv, metadata, final_csv)

        print("")
        print("All ADMET results retrieved.")
        print(f"ADMET molecules prepared: {len(selected)}")
        print("Starting PK Focus extraction for every molecule.")

        page.goto(result_url, wait_until="domcontentloaded",)
        focus_rows, focus_detailed_rows = collect_focus_for_all_molecules(page, output, selected, metadata)
        print("")
        print("==========================================")
        print("DEEP-PK ADMET + PK FOCUS COMPLETED")
        print("==========================================")
        print(f"Molecules: {len(selected)}")
        print(f"PK Focus molecule rows: {focus_rows}")
        print(f"PK Focus detailed rows: {focus_detailed_rows}")
        print(f"Result URL: {result_url}")
        print(f"ADMET CSV: {final_csv}")
        print(f"PK Focus CSV: {output / 'PK_Focus_all_molecules.csv'}")
        print("==========================================")

        browser.close()

if __name__ == "__main__":
    main()