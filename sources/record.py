import numpy as np
from sklearn.cluster import dbscan, OPTICS
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import Angle
import re
from datetime import datetime
from . import db
from astropy.utils.data import compute_hash
import sep
from astropy import wcs
from sdi.sources.reference import reference
import pyvo as vo
from scipy.optimize import curve_fit
import sys

def record(image, path, secid, residual_data, temp, db_session=None):
	"""
	only works with lco shit for now
	"""
	# bkg = sep.background(image.data)
	# recarray = sep.extract(image.data - bkg.back(), bkg.globalrms * 3.0)
	session = db_session
	if session is None:
		session = db.create_session("/seti_data/sdi.williamtest.db")
    
	cat = image["CAT"]
	sci = image["SCI"]
     
	hash_ = compute_hash(path)
	img = session.query(db.Image).filter(db.Image.hash==hash_).first()

	if img is None:
		print(path)
		img = db.Image(image, secid, path)
	
	service = vo.dal.SCSService("https://heasarc.gsfc.nasa.gov/cgi-bin/vo/cone/coneGet.pl?table=m31stars&") 
	w = wcs.WCS(sci.header)
	for element in residual_data:
		trans = db.Transient(data=element, w=w)
		session.add(trans)
		img.transients.append(trans)
		temp.transients.append(trans)
		pixarray = np.array([[element['x'], element['y']]])
		radec = w.wcs_pix2world(pixarray,0)
		result = service.search((radec[0][0], radec[0][1]), 0.001)
		if (result):
			ref = db.Reference(result)
			session.add(ref)
			ref.transients.append(trans)
	appmags = []
	fluxes = []
	source_list = []		
	for source in cat.data:
		# rec = session.query(db.Record).filter(db.Record.ra==r, db.Record.dec==d).first()
		# if rec is None:
			# rec = db.Record(ra=r, dec=d)
	#     session.add(rec)
		s = db.Source(data=source, dtype=cat.data.dtype)
		session.add(s)
		# rec.sources.append(s)
		img.sources.append(s)
		temp.sources.append(s)
		result = reference(source)
		if(result):
			ref = db.Reference(result[0][1])
			session.add(ref)
			ref.sources.append(s)
			appmags.append(result[0][1]['appmag'][0])
			fluxes.append(source[6])
		source_list.append(s)
	coeff, pcov = curve_fit(normalize, fluxes, appmags)
	img.coeff_a = coeff[0]
	img.coeff_b = coeff[1]
	for elem in source_list:
		elem.appmag = coeff[0]*np.log(elem.flux) + coeff[1]
        	
	session.commit()

def _norm(array):
	array -= min(array)
	array *= 1/max(array)
	return array

def flux_sort(source):
	return source[3]

def cluster(db_session=None, algorithm='dbscan', eps=0.001, minpts=4, max_eps=np.inf, xi=0.8, min_flux=0, max_flux=sys.maxsize):
	"""
	Clusters sources using either dbscan or OPTICS from a sqlalchemy session before committing them as records.
	dbscan takes eps and minpts as parameters. OPTICS takes max_eps, minpts, and xi as parameters.
	"""
	session = db_session
	if session is None:
		session = db.create_session()
	irdf = session.query(db.Source.id,db.Source.ra,db.Source.dec,db.Source.flux).all()
	irdf.sort(key=flux_sort)

	# Choose source cutoff based on min/max flux
	min_index = 0
	max_index = len(irdf)
	for i in range(0, len(irdf)):
		if irdf[i][3] > min_flux:
			min_index = i
			break
	for i in range(0, len(irdf)):
		if irdf[-i-1][3] < max_flux:
			max_index = len(irdf)-i
			break
	irdf = np.array(irdf).T

	# do the norming with numpy
	irdf = np.vstack((irdf[0][min_index:max_index],_norm(irdf[1])[min_index:max_index], _norm(irdf[2])[min_index:max_index], _norm(irdf[3])[min_index:max_index]))
	if algorithm == 'dbscan':
		cores, labels = dbscan((irdf[1:]).T, eps, minpts)
	else:
		labels = OPTICS(min_pts, max_eps, xi=xi).fit(irdf[1:].T).labels_
	records = session.query(db.Record).count()
	if records == 0:
		labels += 1 # no -1 label
	else:
		labels += records # Used when running dbscan more than once
	print(irdf[0].shape)
	print(labels.shape)
	for id_, ell in zip(irdf[0], labels):
		print(ell)
		rec = session.query(db.Record).filter(db.Record.label==ell).first()
		if not rec:
			rec = db.Record(label=ell, ra_avg=0, dec_avg=0, flux_avg=0, ra_std=0, dec_std=0, flux_std=0)
		rec.sources.append(session.query(db.Source).get(id_))
		session.add(rec)
	records = session.query(db.Record).all()
	for elem in records:
		sources = elem.sources.all()
		sum_ra = 0
		sum_dec = 0
		sum_flux = 0
		for item in sources:
			sum_ra += item.ra
			sum_dec += item.dec
			sum_flux += item.flux
		elem.ra_avg = sum_ra/len(sources)
		elem.dec_avg = sum_dec/len(sources)
		elem.flux_avg = sum_flux/len(sources)
	session.commit()

def normalize(x, a, b):
	return a*np.log(x) + b

