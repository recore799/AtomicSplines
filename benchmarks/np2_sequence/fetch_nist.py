import urllib.request
import csv
import io

elements = ['C I', 'Si I', 'Ge I']

for el in elements:
    url = f"https://physics.nist.gov/cgi-bin/ASD/energy1.pl?spectrum={el.replace(' ', '+')}&units=1&format=2&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=1&j_out=on&g_out=on&lande_out=on"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req) as response:
            data = response.read().decode('utf-8')
            
            print(f"\n--- NIST Data for {el} ---")
            reader = csv.reader(io.StringIO(data))
            
            headers = []
            for row in reader:
                if row and 'Configuration' in row[0]:
                    headers = row
                    break
                    
            if not headers:
                print("Could not find headers.")
                continue
                
            for row in reader:
                if not row or len(row) < 5: continue
                conf = row[0].strip()
                term = row[1].strip()
                j = row[2].strip()
                level = row[3].strip()
                g = row[4].strip()
                
                if '2s2.2p2' in conf or '3s2.3p2' in conf or '4s2.4p2' in conf or '2s2 2p2' in conf or '3s2 3p2' in conf or '4s2 4p2' in conf or '2p2' in conf or '3p2' in conf or '4p2' in conf:
                    print(f"Conf: {conf:<15} Term: {term:<5} J: {j:<3} Level(cm-1): {level:<15} g: {g}")
                    
    except Exception as e:
        print(f"Error fetching {el}: {e}")
