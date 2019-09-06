
import six
import json

from astropy import units
from xml.dom import minidom
import xml
from astropy.coordinates import SkyCoord

COLUMN_SEPARATOR = "|"


def create_astrores_document():
    implementation = xml.dom.getDOMImplementation()
    doctype = implementation.createDocumentType('ASTRO',
                                                None,
                                                "http://vizier.u-strasbg.fr/xml/astrores.dtd")
    dom = implementation.createDocument("http://vizier.u-strasbg.fr/doc/astrores.htx", "ASTRO", doctype)
    dom.getElementsByTagName("ASTRO")[0].setAttribute("ID", "v0.8")
    dom.getElementsByTagName("ASTRO")[0].setAttribute("xmlns:ASTRO", "http://vizier.u-strasbg.fr/doc/astrores.htx")
    assert isinstance(dom, minidom.Document)
    return dom


class EphemTarget(object):
    nodes = {"NAME": "Ephemeris",
             "TITLE": "Ephemeris for CFHT QSO"}
    fields = {"DATE_UTC": {"attr": {"datatype": "A", "width": "19", "format": "YYYY-MM-DD hh:mm:ss"},
                           "DESCRIPTION": "UTC Date"},
              "RA_J2000": {"attr": {"datatype": "A", "width": "11", "format": "RAh:RAm:RAs", "unit": "h"},
                           "DESCRIPTION": "Right ascension of target"},
              "DEC_J2000": {"attr": {"datatype": "A", "width": "11", "format": "DEd:DEm:DEs", "unit": "deg"},
                            "DESCRIPTION": "Declination of target"}}

    def __init__(self, name, column_separator=COLUMN_SEPARATOR, ephem_format='CFHT ET', runid='16BP06'):
        """
        create an ephmeris target, either with a 'orbfit' object or some mean rate of motion.

        :param name: a string containing the name of the target.
        """

        self.name = str(name).replace(" ", "_")
        self.format = ephem_format
        self.column_separator = column_separator
        self.coordinates = []
        self.runid = runid

    def _init_cfht_api_(self):
        return {'runid': self.runid,
                "pi_login": "mwilson",
                "program_configuration": {"target": []}}

    def _init_cfht_et_file(self):
        self.doc = create_astrores_document()
        doc_root = self.doc.documentElement
        table = self.doc.createElement("TABLE")
        table.setAttribute("ID", "Table")
        doc_root.appendChild(table)

        self.nodes = EphemTarget.nodes
        self.fields = EphemTarget.fields
        nodes = self.nodes
        fields = self.fields
        self.field_names = ["DATE_UTC", "RA_J2000", "DEC_J2000"]

        for (key, value) in six.iteritems(nodes):
            element = self.doc.createElement(key)
            element.appendChild(self.doc.createTextNode(value))
            table.appendChild(element)

        table.getElementsByTagName("TITLE")[0].lastChild.appendData(" target {}".format(self.name))

        table.appendChild(self.doc.createComment("Definition of each field"))

        for fieldName in self.field_names:
            field = self.doc.createElement("FIELD")
            field.setAttribute("name", fieldName)
            for (key, value) in six.iteritems(fields[fieldName]['attr']):
                field.setAttribute(key, value)
            description = self.doc.createElement("DESCRIPTION")
            description.appendChild(self.doc.createTextNode(fields[fieldName]['DESCRIPTION']))
            field.appendChild(description)
            table.appendChild(field)

        table.appendChild(self.doc.createComment("Data table"))

        data = self.doc.createElement("DATA")
        table.appendChild(data)

        header_lines = self._cdata_header(colsep=self.column_separator)

        csv = self.doc.createElement("CSV")
        csv.setAttribute("headlines", str(len(header_lines)))
        csv.setAttribute("colsep", self.column_separator)
        data.appendChild(csv)
        self.cdata = self.doc.createCDATASection("\n" + "\n".join(header_lines) + "\n")
        csv.appendChild(self.cdata)

    @staticmethod
    def _entry(value, width, colsep):
        return "{value:{width}.{width}}{colsep}".format(value=value, width=width, colsep=colsep)

    def _cdata_header(self, colsep="|"):
        """
        Create a header for the CDATA section, as a visual guide.
        """
        fields = self.fields
        header_lines = []
        line = ""
        for fieldName in self.field_names:
            width = int(fields[fieldName]['attr']['width'])
            line += self._entry(fieldName, width, colsep)
        header_lines.append(line)

        line = ""
        for fieldName in self.field_names:
            width = int(fields[fieldName]['attr']['width'])
            line += self._entry(fields[fieldName]['attr']['format'], width=width, colsep=colsep)
        header_lines.append(line)

        line = ""
        for fieldName in self.field_names:
            width = int(fields[fieldName]['attr']['width'])
            (l, m) = divmod(width, 10)
            guide = ""
            for i in range(l):
                guide += "".join(map(str, list(range(10))))
            guide += "".join(map(str, list(range(m))))
            line += self._entry(guide, width=width, colsep=colsep)
        header_lines.append(line)

        line = ""
        for fieldName in self.field_names:
            width = int(fields[fieldName]['attr']['width'])
            guide = "-" * width
            line += self._entry(guide, width=width, colsep=colsep)
        header_lines.append(line)

        return header_lines

    def append(self, coordinate):
        self.coordinates.append(coordinate)

    def _append_cdata(self, coordinate):
        """
        Append an target location to the ephemeris listing.
        """
        fields = self.fields
        sra = coordinate.ra.to_string(units.hour, sep=':', precision=2, pad=True)
        sdec = coordinate.dec.to_string(units.degree, sep=':', precision=1, alwayssign=True)
        coord = SkyCoord(sra + " " + sdec, unit=(units.hour, units.degree))
        sra = coord.ra.to_string(units.hour, sep=":", precision=2, pad=True)
        sdec = coord.dec.to_string(units.degree, sep=":", precision=1, pad=True, alwayssign=True)
        sdate = str(coordinate.obstime.replicate(format('iso')))
        self.cdata.appendData(self._entry(sdate, fields["DATE_UTC"]['attr']['width'], colsep=self.column_separator))
        self.cdata.appendData(self._entry(sra, fields["RA_J2000"]['attr']['width'], colsep=self.column_separator))
        self.cdata.appendData(self._entry(sdec, fields["DEC_J2000"]["attr"]["width"], colsep=self.column_separator))
        self.cdata.appendData("\n")

    def cfht_api_writer(self, f_handle):
        ephemeris_points = []
        for coordinate in self.coordinates:
            epoch_millis = "{:.5f}".format(coordinate.obstime.mjd)
            this_coordinate = {"ra": "{:.4f}".format(float(coordinate.ra.degree)),
                               "dec": "{:.4f}".format(float(coordinate.dec.degree))}
            ephemeris_points.append({"epoch_millis": epoch_millis,
                                     "mag": coordinate.mag,
                                     "coordinate": this_coordinate})
        target = {"identifier": {"client_token": "{}-{}".format(self.runid, self.name)},
                  "name": self.name,
                  "moving_target": {"ephemeris_points": ephemeris_points}}
        json.dump(target, f_handle)
        return

    def gemini_writer(self, f_handle):
        """
        Write out a GEMINI formated OT ephemeris.  This is just a hack of SSD Horizons output.
        """
        f_handle.write(GEMINI_HEADER)
        # Date__(UT)__HR:MN Date_________JDUT     R.A.___(ICRF/J2000.0)___DEC dRA*cosD d(DEC)/dt
        #          1         2         3         4         5         6         7         8         9
        # 123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        # ' 2019-Jan-30 00:00     01 46 56.46 +10 28 54.9 01 47 56.17 +10 34 27.6    3.520

        for coordinate in self.coordinates:
            date = coordinate.obstime.datetime.strftime('%Y-%b-%d %H:%M')[:17]
            f_handle.write(" {:16} {:17.9f}    {:27} {:+8.5f} {:+8.5f}\n".format(date,
                                                                                 coordinate.obstime.jd,
                                                                                 coordinate.to_string('hmsdms',
                                                                                                      sep=' ',
                                                                                                      precision=4,
                                                                                                      pad=True)[:27],
                                                                                 float(0.0),
                                                                                 float(0.0)),
                           )
        f_handle.write(GEMINI_FOOTER)
        return

    def cfht_writer(self, f_handle):
        self._init_cfht_et_file()
        for coordinate in self.coordinates:
            self._append_cdata(coordinate)
        self.doc.writexml(f_handle, indent="  ", addindent="  ", newl='\n')

    def writer(self, f_handle):
        if self.format == 'CFHT ET':
            self.cfht_writer(f_handle)
        elif self.format == "CFHT API":
            self.cfht_api_writer(f_handle)
        elif self.format == 'GEMINI ET':
            self.gemini_writer(f_handle)
        else:
            raise ValueError("{} unknown ET Format".format(self.format))

    def save(self, filename=None):
        if filename is None:
            if self.format == 'CFHT ET':
                filename = "ET_" + self.name + ".xml"
            elif self.format == 'GEMINI ET':
                filename = self.name + ".eph"
            else:
                filename = self.name + ".txt"
        with open(filename, 'w') as f_handle:
            self.writer(f_handle)


GEMINI_HEADER = """*******************************************************************************
 Revised:  Jan 2, 2013                         Titan                       000
                         http://ssd.jpl.nasa.gov/?sat_phys_par
                           http://ssd.jpl.nasa.gov/?sat_elem
 SATELLITE PHYSICAL PROPERTIES:
  Mean Radius (km)       = 2575.5   +-  2.0  Density (g/cm^3) =  1.880 +- 0.004
  Mass (10^22 g)         = 13455.3           Geometric Albedo =  0.2 
  GM (km^3/s^2)          = 8978.13  +-  0.06  V(1,0)          = -1.2 
 SATELLITE ORBITAL DATA:
  Semi-major axis, a (km)= 1221.87 (10^3)  Orbital period     = 15.945421 d
  Eccentricity, e        = 0.0288          Rotational period  = 
  Inclination, i  (deg)  = 0.28
*******************************************************************************
 
 
*******************************************************************************
Ephemeris / WWW_USER Wed Oct 31 17:11:34 2012 Pasadena, USA      / Horizons    
*******************************************************************************
Target body name: Titan (606)                     {source: SAT351}
Center body name: Earth (399)                     {source: DE405}
Center-site name: Mauna Kea
*******************************************************************************
Start time      : A.D. 2012-Nov-01 00:00:00.0000 UT      
Stop  time      : A.D. 2013-Nov-01 00:00:00.0000 UT      
Step-size       : 60 minutes
*******************************************************************************
Target pole/equ : IAU_TITAN                       {East-longitude -}
Target radii    : 2575.0 x 2575.0 x 2575.0 km     {Equator, meridian, pole}    
Center geodetic : 204.527800,19.8261152,4.2078485 {E-lon(deg),Lat(deg),Alt(km)}
Center cylindric: 204.527800,6006.35451,2151.0229 {E-lon(deg),Dxy(km),Dz(km)}
Center pole/equ : High-precision EOP model        {East-longitude +}
Center radii    : 6378.1 x 6378.1 x 6356.8 km     {Equator, meridian, pole}    
Target primary  : Saturn                          {source: DE405}
Interfering body: MOON (Req= 1737.400) km         {source: DE405}
Deflecting body : Sun, EARTH                      {source: DE405}
Deflecting GMs  : 1.3271E+11, 3.9860E+05 km^3/s^2                              
Atmos refraction: NO (AIRLESS)
RA format       : HMS
Time format     : BOTH
EOP file        : eop.121031.p130122                                           
EOP coverage    : DATA-BASED 1962-JAN-20 TO 2012-OCT-31. PREDICTS-> 2013-JAN-21
Units conversion: 1 AU= 149597870.691 km, c= 299792.458 km/s, 1 day= 86400.0 s 
Table cut-offs 1: Elevation (-90.0deg=NO ),Airmass (< 2.000=YES),Daylight (YES)
Table cut-offs 2: Solar Elongation (  0.0,180.0=NO )                           
***************************************************************************************
 Date__(UT)__HR:MN Date_________JDUT     R.A.___(ICRF/J2000.0)___DEC dRA*cosD d(DEC)/dt
***************************************************************************************
$$SOE
"""
GEMINI_FOOTER = """$$EOE
*******************************************************************************
Column meaning:
 
TIME

  Prior to 1962, times are UT1. Dates thereafter are UTC. Any 'b' symbol in
the 1st-column denotes a B.C. date. First-column blank (" ") denotes an A.D.
date. Calendar dates prior to 1582-Oct-15 are in the Julian calendar system.
Later calendar dates are in the Gregorian system.

  The uniform Coordinate Time scale is used internally. Conversion between
CT and the selected non-uniform UT output scale has not been determined for
UTC times after the next July or January 1st.  The last known leap-second
is used over any future interval.

  NOTE: "n.a." in output means quantity "not available" at the print-time.
 
 R.A._(ICRF/J2000.0)_DEC =
   J2000.0 astrometric right ascension and declination of target center.
Corrected for light-time. Units: HMS (HH MM SS.ff) and DMS (DD MM SS.f)
 
 R.A._(a-apparent)__DEC. =
   Airless apparent right ascension and declination of the target center with
respect to the Earth true-equator and the meridian containing the Earth true
equinox of date.  Corrected for light-time, gravitational deflection of light,
stellar aberration, precession & nutation.
   Units: HMS (HH MM SS.ff) and DMS (DD MM SS.f)
 
 Ang-diam =
   The equatorial angular width of the target body full disk, if it were
fully visible to the observer.  Units: ARCSECONDS


 Computations by ...
     Solar System Dynamics Group, Horizons On-Line Ephemeris System
     4800 Oak Grove Drive, Jet Propulsion Laboratory
     Pasadena, CA  91109   USA
     Information: http://ssd.jpl.nasa.gov/
     Connect    : telnet://ssd.jpl.nasa.gov:6775  (via browser)
                  telnet ssd.jpl.nasa.gov 6775    (via command-line)
     Author     : Jon.Giorgini@jpl.nasa.gov

*******************************************************************************
"""
