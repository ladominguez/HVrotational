from HinetPy import Client, win32

# You need a Hi-net account to access their data
client = Client("ladominguez", "love2Smile")

data, ctable = client.get_waveform('0101', '201001010000', 20)
# The request and downloading process will take several minutes
# waiting data request ...
# waiting data downloading ...

win32.extract_sac(data, ctable)


win32.extract_pz(ctable)

