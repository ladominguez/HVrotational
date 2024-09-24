from obspy import read_inventory
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
fdsn_client = Client('IRIS')
t1 = UTCDateTime("2010-09-3T16:30:00.000")
t2 = UTCDateTime("2010-09-3T17:00:00.000")

inv1 = fdsn_client.get_stations(
    network='NZ', station='BFZ', location='10', channel='HHZ',
    starttime=t1, endtime=t2, level='response')

inv_tr130 = read_inventory("reftek151B.resp")
inv_G40T = read_inventory("guralp40T.res")
inv_Sil = read_inventory("sillicon.resp")
inv_Sil[0][0][0].response.response_stages[0].decimation_factor = 1
inv_Sil[0][0][0].response.response_stages[0].decimation_input_sample_rate=100
inv_Sil[0][0][0].plot(0.001, output="DISP", show=True)
