import urllib.request
import re
import os

# Itera a través de los años y los meses de la URL
for year in range(2018, 2024):
    for month in range(1, 13):
        # Construye la URL para el año y el mes actual
        url = f'http://ftp.cptec.inpe.br/goes/goes16/goes16_web/glm_acumulado_nc/{year:04}/{month:02}/'
        
        # Crea la carpeta para el año actual si no existe
        if not os.path.exists(str(year)):
            os.mkdir(str(year))
        
        # Crea la subcarpeta para el mes actual dentro de la carpeta del año si no existe
        month_name = '{:02d}'.format(month)
        month_folder = os.path.join(str(year), month_name)
        if not os.path.exists(month_folder):
            os.mkdir(month_folder)
        
        # Descarga la página HTML y obtiene una lista de los nombres de archivo con extensión .nc
        with urllib.request.urlopen(url) as response:
            html = response.read().decode()
            filenames = re.findall('href="(.*?\.nc)"', html)

        # Descarga cada archivo en la lista en la subcarpeta correspondiente
        for filename in filenames:
            file_url = url + filename
            file_path = os.path.join(month_folder, filename)
            
            if not os.path.exists(file_path):
                print(f'Descargando {file_url}...')
                # urllib.request.urlretrieve(file_url, file_path)

                while True:
                    try:
                        print('Intentando ...')
                        urllib.request.urlretrieve(file_url, file_path)
                        break  # si se descarga exitosamente, rompe el ciclo
                    except urllib.error.URLError:
                        print("Error al descargar el archivo. Intentando nuevamente...")