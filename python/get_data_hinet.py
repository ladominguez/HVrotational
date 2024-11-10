from HinetPy import Client, win32

# You need a Hi-net account to access their data
client = Client("ladominguez", "love2Smile")
client.select_stations("0101", ["N.KWSH", "N.IWNH","N.ASAH","N.NKWH", "N.YHTH"])
data, ctable = client.get_continuous_waveform('0101', '201001010000', 20)
# The request and downloading process will take several minutes
# waiting data request ...
# waiting data downloading ...

win32.extract_sac(data, ctable)


win32.extract_sacpz(ctable)

