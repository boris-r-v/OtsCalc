package ots.calc
import ots.complex.*

typealias Real = Complex
/**
 * Класс пути
 *
 * @param name имя пути нужно только вывода в отладочные сообщения
 * @param mesh сетки для данного пути
 * @param r функция погонного сопротивления рельсов вдоль пути.
 * @param rp функция переходного сопротивления рельсы-земля вдоль пути
 * @param fot подключенные к данному пути отсосы тяговых подстанций
 * @param eps ЭПС находящихся на данном пути
 * @param iRv0 параметр конструктора класса нужен для определения задано ли сопротивления в начале рельсовой линии или его расчитать нужно
 * @param iRvn параметр конструктора класса нужен для определения задано ли значение сопротивления в конце рельсовой линии или его расчитать нужно
 * @param iclU напряжения наводимые в рельсах пути от контакных проводов всех подвесок (пустой на постоянном токе)
 * @param zaz таблица ЗАЗ находящихся на данном пути
 * @param Rtch таблица сосредоточенных точечных сопротивлений в рельсах на данном пути
 * @property Rv0 волновое сопротивление слева Ом
 * @property Rvn волновое сопротивление справа Ом
 * @property m3db трехдиагональнаыя матрица коэффициентов для этого пути
 * @property vectorB массив значений свободных членов
 * @property U значение напряжений  рельс-земля в узлах сетки
 * @property I значения тока в рельсах в узлах сетки
 * @property Ignd значения токов стекающих в землю в узлах метки
 * @property clU напряжения наводимые в рельсах пути от контакных проводов всех подвесок, в узлах сетки
 * @property rlU напряжения наводимые в рельсах от тока из других рельсов, неизвестен на начало расчета
 */
class Track(
    val name: String,
    val mesh: Mesh,
    val r: Array<PV>,
    val rp: Array<PV>,
    val fot: Array<PV>,
    var eps: Array<PV>,
    iRv0 : Real? = null,
    iRvn : Real? = null,
    private val iclU: Array<PV>? = null,
    val zaz: Array<PV> = arrayOf<PV>(),
    val Rtch: Array<PV> = arrayOf<PV>(),
    )
{

    private val rDs = mesh.distributeFunctionOverMesh( r )
    private val rpDs = mesh.distributeFunctionOverMesh( rp )
    internal val Rv0: Real = iRv0 ?: sqrt(rDs.first() * rpDs.first())
    internal val Rvn: Real = iRvn ?: sqrt(rDs.last() * rpDs.last())
    internal var clU: Array<Real> = Array(mesh.size){0.R}
    internal var rlU: Array<Real> = Array(mesh.size){0.R}
    internal val m3db: Array<Array<Real>> = mesh.create3diagMatrixBand(this)
    internal var vectorB: Array<Real> = Array(mesh.size){0.R}
    internal var U: Array<Real> = Array(mesh.size){0.R}
    internal var I: Array<Real> = Array(mesh.size){0.R}
    internal var Ignd: Array<Real> = Array(mesh.size){0.R}
    /**
     * Свойства для хранение исторических данных
     * Массивы U, I после расчета одной мгновенной схемы сохраняются для последующей обработки
     * Сл-но запаситесь оперативой
     */
    internal val histU: Array<ArrayList<Real>> = Array(mesh.size){ arrayListOf() }
    internal val histI: Array<ArrayList<Real>> = Array(mesh.size){ arrayListOf() }
    init {
        mesh.addTrack(this)
    }

    /**
     * Копирует данные из массивов U,I в массивы исторических значения
     */
    fun copy2Hist(){
        U.forEachIndexed { i, value ->
            histU[i].add( value )
        }
        I.forEachIndexed { i, value ->
            histI[i].add( value )
        }
    }
    /**
     * Функция установки массива наведенных напряжений от контакной подвеки на рельсы из исходных данных после расчета матрицы влияний МПС
     */
    fun setClU(){
        //internal var clU: Array<Real> = if (iclU != null ) mesh.distributeFunctionOverMesh( iclU ) else Array(mesh.size()){0.R}
        if ( iclU != null ){
            clU = mesh.distributeFunctionOverMesh( iclU )
        }
    }

    /**
     * Устанавливает новые значения наведенных напряжений на рельсы от контактной сети
     * //FIX ME использвоать в сеттер свойства если возможно ставить из другого типа данных Array<PV>
     */
    fun setClU(up: Array<PV>){
        clU = mesh.distributeFunctionOverMesh( up )
    }
}

/**
 * Заготовка под различные проводимости
 *      свойство пути: private val rlRdp: MutableMap<Track, Array<Real>>
 *      используем так: tr.rlU = tr.rlU + tr2.I * tr2.rlRdp[tr]
 * Вопросы: как задать взаимные междупутные сопротивления в конструкторе
 *      1. Как идентифицировать между какими путями заданы проводимости
 *  Можно проводить идентификацию через имена путей, но тогда они должны быть заданы из исходных данных
 *      2.
 *
 */