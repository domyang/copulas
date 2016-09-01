from datetime import datetime, timedelta

###
# Get daylight saving time (the datetime of that year)
###
# rule: 2nd Sunday in March
def get_daylight_saving_time(year):
    dls_day = datetime.strptime('3/8/'+year, '%m/%d/%Y')
    for i in range(20):
        if dls_day.month == 3 and dls_day.weekday() == 6 and week_of_month(dls_day) == 2:
            break
        dls_day += timedelta(days=1)
    dls_day += timedelta(hours=9) 
    return dls_day
    
###
# Get end of daylight saving (the datetime of that year)
###
# rule: 1st Sunday in November
def get_daylight_saving_end(year):
    dls_end = datetime.strptime('11/1/'+year, '%m/%d/%Y')
    for i in range(20):
        if dls_end.month == 11 and dls_end.weekday() == 6 and week_of_month(dls_end) == 1:
            break
        dls_end += timedelta(days=1)
    dls_end += timedelta(hours=8) 
    return dls_end

###
# Get the week of the month of the given date
###
def week_of_month(date):
    month = date.month
    week = 0
    while date.month == month:
        week += 1
        date -= timedelta(days=7)
    return week
