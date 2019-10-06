import matplotlib.pyplot as plt

def compare_plot( col ):

    """
        This is just a function that will compare the 'I' and 'O'
        dictionaries from a column object. Just call 

            compare_plot( <column object> )

        and it will show you a comparison for a quick check.
    """

    #Define figure
    plt.figure()

    #Ozone
    plt.subplot( 1, 2, 1 )
    plt.semilogx( col.I['lyO'], col.D['lyZ'], 'k-' )
    plt.semilogx( col.O['lyO'], col.D['lyZ'], 'r--' )

    #Temp
    plt.subplot( 1, 2, 2 )
    plt.plot( col.I['lyT'], col.D['lyZ'], 'k-' )
    plt.plot( col.O['lyT'], col.D['lyZ'], 'r--' )

    #Show it
    plt.show()
