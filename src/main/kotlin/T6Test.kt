import ots.calc.*
import ots.calc.Mesh
import ots.calc.Track
import ots.complex.Complex
import ots.complex.R
import ots.complex.complex
import ots.complex.mod
import java.text.DecimalFormat
import java.text.NumberFormat
import java.util.*

private const val empty = Double.NaN

/** Координата ЭПС, км. */
private const val X_EPS = 27.845

/** Координата левого ФОТ, км. */
private const val X_FOT_R = 50.0

/** Координата правого ФОТ. */
private const val X_FOT_L = 0.0

/** Координата справа от правой границы схемы - для функций погонных и переходных сопротивлений. */
private const val X_RIGHT = 100.0

/**
 * Контрольный пример: 1 МПЗ, только основной ход, до шести гл. путей.
 * @param useIclU использовать значения наведенных напряжений.
 * @param useMR использовать значения взаимных сопротвилений между путями.
 * @param addMpss добавить МПС в точки присоединения ФОТ. Будут соединены все пути 1->2, 2->3, 3->4 и т.д.
 * @param trackQty количество главных путей [[1, 6]].
 */
class T6Test(useIclU: Boolean, useMR: Boolean=true, addMpss: Boolean = false, trackQty: Int = 6) {

    init {
        require(trackQty in 1..6)
    }
    // Значения токов ФОТ сделал в соответвии с комсол. Сумма токов ФОТ = току ЭПС
    private val eps = arrayOf(PV(X_EPS, complex(236.0, -185.0)))
    private val fot = arrayOf(PV(X_FOT_L, complex(110.4, -85.679)), PV(X_FOT_R, complex(125.6, -99.321))) // старые значения  complex(144.49, 60.99)    complex(134.08, 50.8)
    private val mesh = Mesh(X_FOT_L-2.0, X_FOT_R+2.0, 0.1)
    private val r = arrayOf(PV(X_RIGHT, complex(0.15, 0.705)))
    private val rp = arrayOf(PV(X_RIGHT, complex(10.0, 0.0)))

    /** Удельные ЭДС индукции в эквивалентных рельсах от токов в подвесках, В/км. */
    private val iclUs = arrayOf(
        // Старые значения закоментировалл
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(34.864, 32.744)), PV(X_FOT_R, complex(-33.598, -30.175)),PV(X_FOT_R+0.01, 0.0.R)),
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(34.237, 31.988)), PV(X_FOT_R, complex(-32.581, -28.836)),PV(X_FOT_R+0.01, 0.0.R)),
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(31.577, 28.647)), PV(X_FOT_R, complex(-29.049, -24.288)),PV(X_FOT_R+0.01, 0.0.R)),
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(30.399, 27.134)), PV(X_FOT_R, complex(-27.542, -22.380)),PV(X_FOT_R+0.01, 0.0.R)),
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(28.498, 24.655)), PV(X_FOT_R, complex(-25.324, -19.612)),PV(X_FOT_R+0.01, 0.0.R)),
        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(27.641, 23.524)), PV(X_FOT_R, complex(-24.487, -18.584)),PV(X_FOT_R+0.01, 0.0.R))

        // Теперь значения наведенки от КС соответствуют комсолу
//        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(-30.080, -29.340)), PV(X_FOT_R, complex(34.768, 33.311)),PV(X_FOT_R+0.01, 0.0.R)),
//        arrayOf(PV(X_FOT_L, 0.0.R),PV(X_EPS, complex(-29.314, -28.335)), PV(X_FOT_R, complex(33.880, 32.190)),PV(X_FOT_R+0.01, 0.0.R))
    )
    private val tracks = Array(trackQty) { i -> makeTrack(i, useIclU) }

    /** Таблица реактивных составляющих взаимных погонных сопротивлений путей, Ом/км. */
    private val mutualReactiveRailResistivities = arrayOf(
        /*                      0      1      2      3      4      5  */
        /* 0 */ doubleArrayOf(empty, 0.341, 0.273, 0.255, 0.229, 0.220),
        /* 1 */ doubleArrayOf(empty, empty, 0.299, 0.273, 0.241, 0.229),
        /* 2 */ doubleArrayOf(empty, empty, empty, 0.341, 0.273, 0.255),
        /* 3 */ doubleArrayOf(empty, empty, empty, empty, 0.299, 0.273),
        /* 4 */ doubleArrayOf(empty, empty, empty, empty, empty, 0.341)
    )

    private val relres = RelativeResist().also {
        val m = mutableMapOf<MRRKey, Array<PV>>()
        for (i in 0..trackQty - 2) {
            for (j in i + 1 until trackQty) {
                val (mrr, pvs) = makeMrrPair(i, j, useMR)
                m[mrr] = pvs
            }
        }
        it.set(MeshRelativeResist(m))
    }
    private val mpss = fot.flatMap { f ->
        tracks.toList().zipWithNext().map { (t1, t2) -> Mps(t1, t2, f.point, f.point, 0.01.R) }
    }
    private val solver = Compute(
        tracks = tracks,
        mpss = if (addMpss) mpss.toTypedArray() else emptyArray(),
        meshes = arrayOf(mesh),
        relres =  relres
    ).apply {
        println(calcOts())
        check(!errorsAndMessages.solverError) {
            "Ошибка в ходе выполнения расчета ОТС: ${errorsAndMessages.messegSolverError}"
        }
    }


    /**
     * Распечатать результаты расчёта.
     * @param precision количество знаков после десятичной запятой.
     */
    fun printResults(precision: Int = 1) {
        fmt.maximumFractionDigits = precision
        printCoordinates(solver.tracks.first().mesh.X, fmt)
        printDataArrays("amps", solver.tracks.map { it.I }, fmt)
        printDataArrays("volts", solver.tracks.map { it.U }, fmt)
        println("Imps= "+solver.getMpsI().mod().contentDeepToString())
    }

    private val fmt = DecimalFormat.getNumberInstance(Locale("en")).apply {
        maximumFractionDigits = 1; isGroupingUsed = false
    }

    private fun printDataArrays(name: String, arrays: List<Array<Complex>>, fmt: NumberFormat) {
        println("$name = [")
        for (a in arrays) {
            print(" [")
            a.forEach { z -> print(fmt.format(z.mod)); print(", ") }
            print("],\n")
        }
        println(']')
    }

    private fun printCoordinates(xx: Array<Double>, fmt: NumberFormat) {
        print("xx = [")
        xx.forEach { x -> print(fmt.format(x)); print(", ") }
        println("]")
    }

    private fun makeTrack(idx: Int, useIclU: Boolean): Track {
        val f = if (idx == 0) fot else emptyArray()
        val e = if (idx == 0) eps else emptyArray()
        val iclU = if (useIclU) iclUs[idx] else null
        return Track((idx + 1).toString(), mesh, r, rp, f, e, iclU = iclU)
    }

    private fun makeMrrPair(t1: Int, t2: Int, useMR: Boolean): Pair<MRRKey, Array<PV>> {
        val ra = if (useMR) 0.05 else 0.0
        val rr = if (useMR) mutualReactiveRailResistivities[t1][t2] else 0.0
        return Pair(
            MRRKey(tracks[t1], tracks[t2]),
            arrayOf(PV(X_RIGHT, complex(ra, rr)))
        )
    }
}

fun main() {
    val test = T6Test(true, true,true,trackQty = 6) //расчет с двумя путями, МПС и взаимная индуктивность есть
    test.printResults()
}
