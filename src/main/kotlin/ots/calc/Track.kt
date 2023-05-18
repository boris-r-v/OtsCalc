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
 * @param fot таблица подключенных к данному пути ФОТ
 * @param esp таблица ЭПС находящихся на данном пути
 * @param uRv0 параметр конструктора класса нужен для определения задано ли сопротивления в начале рельсовой линии или его расчитать нужно
 * @param uRvn параметр конструктора класса нужен для определения задано ли значение сопротивления в конце рельсовой линии или его расчитать нужно
 * @param zaz таблица ЗАЗ находящихся на данном пути
 * @param Rtch таблица сосредоточенных точечных сопротивлений в рельсах на данном пути
 * @property Rv0 волновое сопротивление слева Ом
 * @property Rvn волновое сопротивление справа Ом
 * @property m3db трехдиагональнаыя матрица коэффициентов для этого пути
 * @property vectorB - массив хначений свободных членов
 * @property U - значение напряжений  рельс-земля в узлах сетки
 * @property I - значения тока в рельсах в узлах сетки
 * @property Ignd - значения токов стекающих в хемлю в узлах метки
 */
class Track(
    val name: String,
    val mesh: Mesh,
    val r: Array<PV>,
    val rp: Array<PV>,
    val fot: Array<PV>,
    var eps: Array<PV>,
    uRv0 : Real? = null,
    uRvn : Real? = null,
    val zaz: Array<PV> = arrayOf<PV>(),
    val Rtch: Array<PV> = arrayOf<PV>(),

    )
{
    private val rDs = mesh.distributeFunctionOverMesh( r )
    private val rpDs = mesh.distributeFunctionOverMesh( rp )
    internal val Rv0: Real = uRv0 ?: sqrt(rDs.first() * rpDs.first())
    internal val Rvn: Real = uRvn ?: sqrt(rDs.last() * rpDs.last())
    internal val m3db: Array<Array<Real>> = mesh.create3diagMatrixBand(this)
    internal var vectorB: Array<Real> = Array(mesh.size()){0.R}
    internal var U: Array<Real> = Array(mesh.size()){0.R}
    internal var I: Array<Real> = Array(mesh.size()){0.R}
    internal var Ignd: Array<Real> = Array(mesh.size()){0.R}

}