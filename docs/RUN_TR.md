# R kodlarını çalıştırma

RStudio'da çalışma dizini olarak `00_RUN.R` ve `RUN_BMA_REVISION.R`
dosyalarının bulunduğu depo kökünü seç.

~~~r
source("INSTALL_PACKAGES.R")
source("00_RUN.R")
main("audit")
main("smoke")
~~~

İlk komut eksik paketleri kurar. `audit` sayısal denetimleri çalıştırır.
`smoke` kısa benzetimi çalıştırıp çıktıları `results/smoke/` içine kaydeder.
Bu iki adım piyasa verisi gerektirmez.

Tam çalışma için sekiz varlığın fiyat dosyalarını `data/raw/` klasörüne koy.
Dosya adları, sütunlar ve tarih aralıkları README'de yer alıyor.

~~~r
source("RUN_BMA_REVISION.R")
bma_revision("full")
~~~

Tam çalışma `runs/` altında bir koşu dizini açar. Kesinti olursa
aynı R ve paket sürümleriyle koşu dizinini vererek devam edebilirsin.

~~~r
bma_revision("resume", "runs/run_example")
~~~

`run_example` yerine kendi koşu klasörünün adını yaz. Fiyat ve tahmin
dosyaları bulunan tamamlanmış bir koşuyu, kendi dizininden
`Rscript RUN_BMA_REVISION.R verify` komutuyla denetleyebilirsin. Arşivlenmiş
koşular için arşivle birlikte verilen denetim betiklerini kullan. Yeni fiyat
indirmesi veya farklı bir optimizasyon sürümü tahminleri değiştirebilir.
