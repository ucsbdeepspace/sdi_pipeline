from astropy.io import fits
import click
from astroplan import Observer
from astropy import units as u
from . import _cli as cli
from astropy.time import Time,TimezoneInfo
from astropy.units import deg, m
from astropy.coordinates import EarthLocation
from datetime import datetime


def middle_night(hduls, read_ext='SCI', write_ext=-1):
    if write_ext == -1:
        write_ext = read_ext
    for hdul in hduls:
        hdr = hdul['SCI'].header['DATE-OBS']
        halfexpos = float(hdul['SCI'].header['EXPTIME'])/2
        
        year = float(hdr[0:4])
        month = float(hdr[5:7])
        day = float(hdr[8:10])

        obs_hour = hdr[11:13]
        obs_minute = hdr[14:16]
        obs_second = float(hdr[17:])+halfexpos

        obs = Observer(latitude=hduls[0]['SCI'].header['LATITUDE']*deg, longitude=hduls[0]['SCI'].header['LONGITUD']*deg, elevation=hduls[0]['SCI'].header['HEIGHT']*m, timezone='UTC')

        solar_midnight=str(obs.midnight(Time(hdr), which='next').to_datetime())

        hour = int(obs_hour) - int(solar_midnight[11:13])
        minute = int(obs_minute) - int(solar_midnight[14:16])
        second = obs_second - float(solar_midnight[17:])
        if minute < int(0):
            minute += int(60)
            hour -= int(1)
        if second < float(0):
            second += float(60)
            minute -= int(1)
        
        a = float((14-month)/12)
        y = float(year + 4800 - a)
        m = float(month + (12*a) - 3)
        JDN = float(day + (((153*m)+2)/5) + (365*y) + ((97*y)/400) - 32045)
        JD  = float(JDN + (((hour - 12)/24) + (minute / 1440) + (second / 86400)))
        MJD = float(JD - 2400000.5)
        hdul[write_ext].header.set('MIDNT', MJD, comment='Difference between DATE-OBS and solar midnight')
        yield hdul

@cli.cli.command("middle_night")
@click.option("-r", "--read_ext", default='SCI',
              help="An index number or ext name that identifies the data in"
              "input hduls that you want source extraction for. For LCO, this "
              "is 0 or SCI.")
@click.option("-w", "--write_ext", default=-1,
              help="An extension name and extension version that will identify"
              "the HDUL that the resulting middle night time gets written to. Default"
              "is -1")
@cli.operator

def middle_night_cmd(hduls, read_ext, write_ext):
    return middle_night(hduls, read_ext, write_ext)


