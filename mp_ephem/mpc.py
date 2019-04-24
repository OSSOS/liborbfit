"""Interact with the minor planet center"""
import requests
import re
import logging
from . import EphemerisReader
import io

MPCSERVER = "https://www.minorplanetcenter.net"

SHOW_OBJECT_ENDPOINT = "/db_search/show_object"


def download_mpc_database(target):
    """
    Retrieve the astrometric measurements available for the given object.
    :param target:
    :return:
    """
    params  = {"object_id": target}
    object_info_page = requests.get("{}{}".format(MPCSERVER, SHOW_OBJECT_ENDPOINT), params=params)

    #   <a href="../tmp/29981.txt" target="_blank">download</a>
    download_match = re.search(".*href=\"(../)?(?P<filename>[^\"]*)\".*download.*", object_info_page.content)
    retrieval_url = "{}/{}".format(MPCSERVER, download_match.group('filename'))
    try:
        return requests.get(retrieval_url).content
    except Exception as ex:
        logging.error(ex)
        return None


def get_mpc_observations(target):

    fobj = io.StringIO(download_mpc_database(target))
    return EphemerisReader().read(fobj)

