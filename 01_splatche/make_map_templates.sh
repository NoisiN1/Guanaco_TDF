#!/bin/bash
cd /home/apena/Splatche/datasetC

# Copiar plantillas base a directorios de mapas
rm -rf mapas_splatche/scen1_map mapas_splatche/scen2_map
cp -r Splatche_scen1_deme_TDF_C mapas_splatche/scen1_map
cp -r Splatche_scen2_hum_TDF_C  mapas_splatche/scen2_map

# Activar outputs espaciales ASCII (density + occupation + migration)
for scen in scen1_map scen2_map; do
  sed -i 's/GenerateOutputMigrationASCII=0/GenerateOutputMigrationASCII=1/' \
    mapas_splatche/$scen/settings.txt
  sed -i 's/GenerateOutputMDensityASCII=0/GenerateOutputMDensityASCII=1/' \
    mapas_splatche/$scen/settings.txt
  sed -i 's/GenerateOutputOccupationASCII=0/GenerateOutputOccupationASCII=1/' \
    mapas_splatche/$scen/settings.txt
done

# Sustituir K y m con los valores del modo posterior:
#   Scen1 (Refugio): K=36, m=0.056
#   Scen2 (Puente): usaremos la mediana del prior como referencia neutra,
#                   ya que scen2 no fue seleccionado; K=22150, m=0.056
sed -i "s/PARAM/36/g" mapas_splatche/scen1_map/datasets/veg2K_pop1_time_3.txt
sed -i "s/MigrationRate=0.01/MigrationRate=0.056/" mapas_splatche/scen1_map/settings.txt

sed -i "s/PARAM/22150/g" mapas_splatche/scen2_map/datasets/veg2K_pop1_time_5.txt
sed -i "s/PARAM/22150/g" mapas_splatche/scen2_map/datasets/veg2K_pop1_time_6.txt
sed -i "s/MigrationRate=0.01/MigrationRate=0.056/" mapas_splatche/scen2_map/settings.txt

# Verificar
echo "--- Scen1 flags ASCII ---"
grep "ASCII=" mapas_splatche/scen1_map/settings.txt
echo "--- Scen1 MigrationRate y K ---"
grep "MigrationRate=" mapas_splatche/scen1_map/settings.txt
grep -c "PARAM" mapas_splatche/scen1_map/datasets/veg2K_pop1_time_3.txt

echo "--- Scen2 flags ASCII ---"
grep "ASCII=" mapas_splatche/scen2_map/settings.txt
