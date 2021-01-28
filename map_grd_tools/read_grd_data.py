# standard libs
import bz2
import datetime
# SD libs
import pydarnio
import pandas
import numpy


class ReadGrdData(object):
    """ 
        Class to Grid level data
        
        Parameters
        ----------
        start_date : datetime
            the start time as a datetime
        end_date : datetime
            the end time as a datetime
        fileType : Optional[str]
            the filetype to read, "grdex","grid2";
    """
    
    def __init__(self, start_date, end_date,\
                 file_type="grid2",\
                 hemi=["north", "south"],\
                convert_to_ascii=True,
                out_dir="/tmp/"):
        """ 
        Initialize all parameters 
        """
        self.start_date = start_date
        self.end_date = end_date
        self.file_type = file_type
        self.convert_to_ascii = convert_to_ascii
        self.out_dir = out_dir
        if isinstance(hemi, str):
            hemi = [hemi]
        self.hemi = hemi
        
        _ndays = (self.start_date - self.end_date).days + 2
        self.file_list = []
        while start_date <= end_date:
            for _hm in self.hemi:
                _filestr = "/sd-data/{year}/{file_type}/{hemi}/{date}.{hemi}.{file_type}.bz2"
                _filestr = _filestr.format(year=start_date.year,\
                                hemi=_hm,file_type=self.file_type,\
                                date=start_date.strftime("%Y%m%d"))
                self.file_list.append(_filestr)
            start_date += datetime.timedelta(days=1)
        
    def save_to_ascii(self):
        """
        Read the data from files
        """
        for _file in self.file_list:
            print("working through-->", _file)
            _fend = _file.split("/")[-1]
            _fname = ".".join(_fend.split(".")[:-1]) + ".txt"
            out_file_name = self.out_dir + "/" + _fname
            outf = open(out_file_name, "w")
            with bz2.open(_file) as fp: 
                  data_stream = fp.read()
            reader = pydarnio.SDarnRead(data_stream, True)
            rec_list = reader.read_grid()
            for _recs in rec_list:
                _stime = datetime.datetime( _recs["start.year"], _recs["start.month"],\
                                           _recs["start.day"], _recs["start.hour"],\
                                           _recs["start.minute"], int(_recs["start.second"]) )
                
                _etime = datetime.datetime( _recs["end.year"], _recs["end.month"],\
                                           _recs["end.day"], _recs["end.hour"],\
                                           _recs["end.minute"], int(_recs["end.second"]) )
                
                
                outf.write("start date = " + _stime.strftime("%Y-%m-%d  ") + "\n")
                outf.write("start time = " + _stime.strftime("%H:%M:%S  ") + "\n")
                outf.write("end date = " + _etime.strftime("%Y-%m-%d  ") + "\n")
                outf.write("end time = " + _etime.strftime("%H:%M:%S  ") + "\n")
                outf.write("***********************Data Summary*********************** \n")
                               
                # Write the column header
                outf.write("{0:>5s} {1:>5s} {2:>5s} {3:>10s} {4:>5s} {5:>5s} "
                        "{6:>5s} {7:>9s} {8:>10s} {9:>9s} {10:>5s} {11:>5s} "
                        "{12:>5s} {13:>5s} {14:>5s} {15:>5s} {16:>5s} {17:>5s}\n".
                        format("stid", "channel", "nvec", "freq", "major.revision", "minor.revision",
                               "program.id", "noise.mean", "noise.sd", "gsct",
                               "v.min", "v.max", "p.min", "p.max", "w.min", "w.max", "ve.min", "ve.max"))
                
                for _ind, _st in enumerate(_recs['stid']):
                    outf.write("{0:>5s} {1:>5s} {2:>5s} {3:>10s} {4:>5s} {5:>5s} "
                        "{6:>5s} {7:>9s} {8:>10s} {9:>9s} {10:>5s} {11:>5s} "
                        "{12:>5s} {13:>5s} {14:>5s} {15:>5s} {16:>5s} {17:>5s}\n".
                            format(str(_st), str(_recs["channel"][_ind]), str(_recs["nvec"][_ind]),\
                           str(_recs["freq"][_ind]),str(_recs["major.revision"][_ind]),\
                           str(_recs["minor.revision"][_ind]), str(_recs["program.id"][_ind]),\
                           str(_recs["noise.mean"][_ind]), str(_recs["noise.sd"][_ind]),\
                           str(_recs["gsct"][_ind]), str(_recs["v.min"][_ind]),\
                           str(_recs["v.max"][_ind]), str(_recs["p.min"][_ind]),\
                           str(_recs["p.max"][_ind]), str(_recs["w.min"][_ind]),\
                           str(_recs["w.max"][_ind]), str(_recs["ve.min"][_ind]),\
                           str(_recs["ve.max"][_ind])))
                    outf.write("\n")
                outf.write("***********************Data Summary*********************** \n")
                
                outf.write("***********************Data Records*********************** \n")
                # Write the column header
                outf.write("{0:>5s} {1:>10s} {2:>11s} {3:>5s} {4:>5s} {5:>5s} "
                                "{6:>11s} {7:>11s} {8:>10s} {9:>10s} {10:>10s} {11:>10s}\n".
                        format("mlat", "mlon", "kvect",\
                               "stid", "channel", "index",\
                               "vel.median", "vel.sd",\
                               "pwr.median", "pwr.sd",\
                               "wdt.median", "wdt.sd"))
                if not 'vector.mlat' in _recs:
                    outf.write("No records \n")
                    print("No records-->", _stime)
                else:
                    for _ind, _mlat in enumerate(_recs['vector.mlat']):
                            outf.write("{0:>5s} {1:>10s} {2:>11s} {3:>5s} {4:>5s} {5:>5s} "
                                "{6:>11s} {7:>11s} {8:>10s} {9:>10s} {10:>10s} {11:>10s}\n".
                                    format(str(_mlat), str(_recs["vector.mlon"][_ind]), str(_recs["vector.kvect"][_ind]),\
                                   str(_recs["vector.stid"][_ind]), str(_recs["vector.channel"][_ind]),\
                                   str(_recs["vector.index"][_ind]),str(_recs["vector.vel.median"][_ind]),\
                                   str(_recs["vector.vel.sd"][_ind]), str(_recs["vector.pwr.median"][_ind]),\
                                   str(_recs["vector.pwr.sd"][_ind]), str(_recs["vector.wdt.median"][_ind]),\
                                   str(_recs["vector.wdt.sd"][_ind])))
                            outf.write("\n")

                outf.write("***********************Data Records*********************** \n")
            outf.close()
            print("saved text to-->", out_file_name)
            
                
            
        
if __name__ == "__main__":
    start_date = datetime.datetime(2017,5,26)
    end_date = datetime.datetime(2017,5,31)
    hemi = ["north", "south"]
    grd_obj = ReadGrdData(start_date, end_date, hemi=hemi, out_dir="./")
    grd_obj.save_to_ascii()