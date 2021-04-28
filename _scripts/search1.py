import db
from glob import glob
from astropy.io import fits
import numpy as np
from vis import image
from loaddb import load_im
import sys
from sdi.sources.record import record, cluster

ims = sys.argv[1:]
#fn = glob("smalldata/*.fz")
#hduls = [fits.open(f) for f in fn]
#sci = [h['SCI'] for h in hduls]
#cat = [h['CAT'] for h in hduls]

s = db.create_session('/home/pkotta/GTAnd_PTF_small.db')
time_id = {}


def id_func(source):
    return source.image_id

def time_func(image):
    return image.time

def get_rec(session, record):
    return(session.query(db.Record).get(record).sources.all())

def test_func(source):
    return source.image.time

def rec_search(ra, dec, session = s):
    min_dif = (session.query(db.Record).get(2).sources.all()[0].ra_avg - ra)**2 + (session.query(db.Record).get(2).sources.all()[0].dec_avg - dec)**2
    rec = session.query(db.Record).get(2)
    for record in session.query(db.Record).all()[2:]:
        dif = (record.sources.all()[0].ra_avg - ra)**2 + (record.sources.all()[0].dec_avg-dec)**2
        if dif < min_dif:
            min_dif = dif
            rec = record
    return rec
        
#Functions for counting gaps

def gap_func(sources):
    sources = sources
    gaps = 0
    for i in range(0, len(sources)-1):
        if sources[i].image_id + 2 == sources[i+1].image_id:
            gaps += 1
    return gaps

def full_gaps(session, just_counts = True):
    gap_list = []
    records = s.query(db.Record).count()
    for i in range(2, records + 1):
        gap_list.append(gap_func(get_rec(session, i)))
    if just_counts == True:
        return sum(gap_list)
    else:
        return gap_list, sum(gap_list)

def gap_im_count(session):
    gap_count = []
    for i in range(2, records+1):
        images = []
        time_ids = []
        sources = session.query(db.Record).get(i).sources.all()
        for source in sources:
            images.append(source.image)
        images.sort(key=time_func)
        for image in images:
            time_ids.append(time_id[str(image.id)])
        for j in range(0, len(time_ids)-1):
            if time_ids[j] + 2 == time_ids[j+1]:
                gap_count.append(images[j].id)
    return gap_count

def gap_im_source(session, im_id):
    source_count = []
    for i in range(2, records + 1):
        time_ids = []
        sources = session.query(db.Record).get(i).sources.all()
        sources.sort(key=test_func)
        for source in sources:
            time_ids.append(time_id[str(source.image.id)])
        for j in range(0, len(time_ids)-1):
            if time_ids[j] + 2 == time_ids[j+1] and sources[j].image_id == im_id:
                source_count.append((sources[j].x, sources[j].y))
                break
    source_count = np.array(source_count)
    return source_count

#Functions for counting collisions
def col_func(sources):
    sources.sort(key = id_func)
    collisions = 0
    for i in range(0, len(sources)-1):
        if sources[i].image_id == sources[i+1].image_id and sources[i].x != sources[i+1].x:
            collisions += 1
    return collisions

def full_col(session, just_counts = True):
    col_list = []
    records = session.query(db.Record).count()
    for i in range(2, records + 1):
        col_list.append(col_func(session.query(db.Record).get(i).sources.all()))
    if just_counts == True:
        return sum(col_list)
    else:
        return col_list, sum(col_list)

def col_im_count(session):
    col_count = []
    for i in range(2, records + 1):
        sources = get_rec(session, i)
        sources.sort(key = id_func)
        for j in range(0, len(sources)-1):
            if sources[j].x == sources[j+1].x:
                col_count.append(sources[j].image_id)
    return col_count

def multi_source_count(session):
    col_count = []
    for i in range(2, records + 1):
        sources = get_rec(session, i)
        sources.sort(key = id_func)
        for j in range(0, len(sources)-1):
            if sources[j].image_id == sources[j+1].image_id and sources[j].x != sources[j+1].x:
                col_count.append(sources[j].image_id)
    return col_count, sum(col_count)

def col_im_sources(session, im_id):
    source_count = []
    source_rec = []
    for i in range(2, records + 1):
        sources = get_rec(session, i)
        sources.sort(key = id_func)
        for j in range(0, len(sources)-1):
            if sources[j].image_id == sources[j+1].image_id and sources[j].image_id == im_id:
                source_count.append((sources[j].x,sources[j].y))
                source_rec.append((sources[j], sources[j+1]))
                break
    source_count = np.array(source_count)
    return source_count, source_rec
"""
#Functions for counting gaps and collisions
def count(session, record, gap_size = 1, show_rec = False):
    sources = get_rec(session, record)
    sources.sort(key = func)
    gaps = 0
    collisions = 0
    for i in range(0, len(sources)-1):
        if sources[i].image_id == sources[i+1].image_id:
            collisions += 1
        elif sources[i+1].image_id >= sources[i].image_id + 2 and sources[i+1].image_id < sources[i].image_id + 2 + gap_size:
            gaps += 1
    if show_rec == True:
        return (gaps, collisions, record)
    else:
        return (gaps, collisions)

#Returns gaps and collisions for all records as a list of tuples, record 1 is skipped
def full_count(session, gap_size=1, show_rec = False):
    records = session.query(db.Record).count()
    counts = []
    for i in range(2, records + 1):
        counts.append(count(session, i, gap_size, show_rec))
    return counts

#First entry/tuple in full count functions corrsepond to Record:2
#e.g., full_counts[i] => Record:i+2
"""

def ctn(session = s):
    noise = 0
    for source in session.query(db.Source).all():
        if source.record.id == 1:
            noise += 1
    return noise/session.query(db.Source).count()

def super_func(session=s):
    print('cols',full_col(session))
    print('gaps',full_gaps(session))
    print("Records:", session.query(db.Record).count())
    print('cluster to noise', ctn(session))

def noise_points(im_id, session = s):
    points = []
    for source in session.query(db.Image).get(im_id).sources.all():
        if source.record.id == 1:
            points.append((source.x, source.y))
    return points

#Clear all entries
def clear_all(session=s):
    for image in session.query(db.Image).all():
        session.delete(image)
    for record in session.query(db.Record).all():
        session.delete(record)
    for source in session.query(db.Source).all():
        session.delete(source)
    session.commit()

def clear_rec(session=s):
    for record in session.query(db.Record).all():
        session.delete(record)

a = []
def do_thing():
    for i in range(4, 11):
        clear_rec(s)
        cluster(s, min_dist = 0.0018, min_points = i)
        a.append([full_col(s), full_gaps(s), s.query(db.Record).count(), i])
    print(a)
