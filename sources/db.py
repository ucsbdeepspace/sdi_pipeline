import sqlalchemy
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column, Integer, DateTime, String, Text, Float, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship, backref
from astropy.utils.data import compute_hash
from astropy.coordinates import Angle
from datetime import datetime
import re
import numpy as np

Base = declarative_base()


class Record(Base):
    __tablename__ = "Record"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    label = Column(Integer, unique=True)
    ra_avg = Column(Float)
    dec_avg = Column(Float)
    flux_avg = Column(Float)
    ra_std = Column(Float)
    dec_std = Column(Float)
    flux_std = Column(Float)
    sources = relationship("Source", backref="record", lazy="dynamic", foreign_keys="Source.record_id")

    def __repr__(self):
        return "<Record {} RA:{} DEC:{} FLUX:{}>".format(self.id, self.ra_avg, self.dec_avg, self.flux_avg)


class Source(Base):
    __tablename__ = "Source"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    image_id = Column(Integer, ForeignKey('Image.id'))
    time = Column(DateTime, ForeignKey('Image.time'))
    record_id = Column(Integer, ForeignKey('Record.id'))
    reference_id = Column(Integer, ForeignKey('Reference.id'))
    template_id = Column(Integer, ForeignKey('Template.id'))
    x = Column(Float)
    y = Column(Float)
    xwin = Column(Float)
    ywin = Column(Float)
    xpeak = Column(Float)
    ypeak = Column(Float)
    flux= Column(Float)
    fluxerr = Column(Float)
    peak = Column(Float)
    fluxaper1 = Column(Float)
    fluxerr1= Column(Float)
    fluxaper2 = Column(Float)
    fluxerr2 = Column(Float)
    fluxaper3 = Column(Float)
    fluxerr3 = Column(Float)
    fluxaper4 = Column(Float)
    fluxerr4 = Column(Float)
    fluxaper5= Column(Float)
    fluxerr5 = Column(Float)
    fluxaper6 = Column(Float)
    fluxerr6 = Column(Float)
    background = Column(Float)
    fwhm = Column(Float)
    a = Column(Float)
    b = Column(Float)
    theta= Column(Float)
    kronrad = Column(Float)
    ellipticity = Column(Float)
    fluxrad = Column(Float)
    fluxrad50 = Column(Float)
    fluxrad75= Column(Float)
    x2 = Column(Float)
    y2 = Column(Float)
    xy= Column(Float)
    flag = Column(Float)
    ra = Column(Float)
    dec = Column(Float)
    appmag = Column(Float)

    def __init__(self, data, dtype=None, kwargs={}):
        """
        Initializes the Source using a line in a numpy recarray
        """
        super(Source, self).__init__(**kwargs)

        if dtype is None:
            dtype = data.dtype
        for idx, n in enumerate(dtype.names):
            if n.lower() in Source.__dict__:
                self.__dict__[n.lower()] = data[idx]

 
    def __repr__(self):
        return "<Source {} Image:{} Record:{}>".format(self.id, self.image_id, self.record_id)

class Transient(Base):
    __tablename__ = "Transient"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    image_id = Column(Integer, ForeignKey('Image.id'))
    reference_id = Column(Integer, ForeignKey('Reference.id'))
    template_id = Column(Integer, ForeignKey('Template.id'))
    thresh = Column(Float)
    ra = Column(Float)
    dec = Column(Float)
    def __init__(self, data, w):
        pixarray = np.array([[data['x'],data['y']]])
        radec = w.wcs_pix2world(pixarray,0)
        self.ra = radec[0][0]
        self.dec = radec[0][1]
        self.thresh = data['thresh']
    def __repr__(self):
        return "<Transient {} RA:{} DEC:{}>".format(self.id, self.ra, self.dec)


class Template(Base):
    __tablename__ = "Template"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    path = Column(Text(255), unique=True)
    section_id = Column(Text(255))
    sources = relationship("Source", backref="template", lazy="dynamic", foreign_keys="Source.template_id")
    transients = relationship("Transient", backref="template", lazy="dynamic", foreign_keys="Transient.template_id") 
    def __init__(self, path, secid):
	
        self.path = path
        self.section_id = secid
	 
    def __repr__(self):
        return "<Template {} Section:{} Path:{}>".format(self.id, self.section_id, self.path)


class Reference(Base):
    __tablename__ = "Reference"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    ra = Column(Float)
    dec = Column(Float)
    name = Column(Text(255))
    lii = Column(Float)
    bii = Column(Float)
    appmag = Column(Float)
    appmag_error = Column(Float)
    type_appmag = Column(Text(255))
    orig_file = Column(Text(255)) 
    sources = relationship("Source", backref="reference", lazy="dynamic", foreign_keys="Source.reference_id")
    transients = relationship("Transient", backref="reference", lazy="dynamic", foreign_keys="Transient.reference_id")
    def __init__(self, ref):
        self.ra = ref['ra'][0]
        self.dec = ref['dec'][0]
        self.name = ref['name'][0]
        self.lii = ref['lii'][0]
        self.bii = ref['bii'][0]
        self.appmag = ref['appmag'][0]
        self.appmag_error = ref['appmag_error'][0]
        self.orig_file = ref['orig_file'][0]
    def __repr__(self):
        return "<Reference {} RA:{} Dec:{}>".format(self.id, self.ra, self.dec)

class Image(Base):
    __tablename__ = "Image"
    id = Column(Integer, primary_key=True, unique=True, nullable=False)
    path = Column(Text(255), unique=True)
    time = Column(DateTime, unique=True)
    ra = Column(Float)
    ra_max = Column(Float)
    ra_min = Column(Float)
    dec = Column(Float)
    coeff_a = Column(Float)
    coeff_b = Column(Float)
    sources = relationship("Source", backref="image", lazy="dynamic", foreign_keys="Source.image_id")
    section_id = Column(Text(255))
    transients = relationship("Transient", backref="Image", lazy="dynamic", foreign_keys="Transient.image_id")
    hash = Column(Text(32), unique=True, index=True)

    def __init__(self, image_or_path, secid, path=None):
        im_path = ""
        img = image_or_path
        if path is None:
            im_path = image_or_path
            img = fits.open(image_or_path)
        else:
            im_path = path

        cat = img["CAT"]
        sci = img["SCI"]
        hash_ = compute_hash(im_path)
        datestr = re.search(r"\d{4}-\d{2}-\d{2}", sci.header["DATE"]).group()
        timestr = re.search(r"\d{2}:\d{2}:\d{2}\.?\d+", sci.header["UTSTART"]).group()
        dt = datetime.strptime(" ".join([datestr, timestr]), "%Y-%m-%d %H:%M:%S.%f")

        self.path = path
        self.time = dt
        self.hash = hash_
        self.ra = Angle(sci.header["RA"], unit="hourangle").deg
        self.dec = Angle(sci.header["DEC"], unit="degree").deg
        self.section_id = secid
    def __repr__(self):
        return "<Image {} Time:{} Path:{}>".format(self.id, self.time, self.path)


def create_session(db_path="/seti_data/sdi.db"):
    engine = create_engine(r"sqlite:///{}".format(db_path), echo=False) #connect to database
    Base.metadata.create_all(engine) #Lets create the actual sqlite database and schema!
    return sessionmaker(bind=engine)()
