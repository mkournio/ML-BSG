import lightkurve as lk

#q = lk.search_targetpixelfile('TIC 90256729',mission="TESS")
q=lk._search_products(['TIC 90256729','TIC 385241538'],filetype="ffi",mission="TESS")
print(q)
q.download_all(download_dir='./')
