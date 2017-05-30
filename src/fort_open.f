      subroutine fortopen(FNAME)
      CHARACTER*25 FNAME

      open(unit=7,file=FNAME,status="unknown")

      end


      subroutine fortopeninput(NAME)
      CHARACTER*25 NAME
      open(unit=1,file=NAME,status="old")
      return
      end




