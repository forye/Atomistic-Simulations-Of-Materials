__author__ = 'Idan'

import datetime
import configuration as conf

def icheck(cond, message,stop_if_fails = True):
    '''
    :param cond: the test itself. returns bool where True == FAILED
    :param message: describes exactly what is the test was
    :param stop_if_fails: if fail (cond == true), throw an exceptions
    '''
    if cond:
        print "->" + message
        open(conf.DEBUG, 'a').write("\n       ---  log from "+ str(datetime.datetime.now().strftime("%H:%M:%S.%f")) +
                                    '---:\n ' + message + '\n')
        if stop_if_fails:
            raise Exception(message)