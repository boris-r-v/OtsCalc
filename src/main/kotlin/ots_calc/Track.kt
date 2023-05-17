package ots_calc
import org.kotlinmath.*
import java.util.Vector

typealias Real = Complex
/**
 * Класс пути
 *
 * @param r функция погонного сопротивления рельсов вдоль пути.
 * @param rp функция переходного сопротивления рельсы-земля вдоль пути
 * @param fot таблица подключенных к данному пути ФОТ
 * @param esp таблица ЭПС находящихся на данном пути
 * @param zaz таблица ЗАЗ находящихся на данном пути
 * @param Rtch таблица сосредоточенных точечных сопротивлений в рельсах на данном пути
 * @param numMesh номер сетки для данного пути
 * @param Rv0 волновое сопротивление слева Ом
 * @param Rvn волновое сопротивление справа Ом
 * @property m3db трехдиагональнаыя матрица коэффициентов для этого пути
 * @property vectorB - массив хначений свободных членов
 * @property U - значение напряжений  рельс-земля в узлах сетки
 * @property I - значения тока в рельсах в узлах сетки
 * @property Ignd - значения токов стекающих в хемлю в узлах метки
 */
class Track(
    val name: String,
    val r: Array<PV>,
    val rp: Array<PV>,
    val fot: Array<PV>,
    val eps: Array<PV>,
    val zaz: Array<PV>,
    val Rtch: Array<PV>,
    val mesh: Mesh,
    var Rv0: Real = -1.R,
    var Rvn: Real = -1.R,
)
{
    internal var m3db: Array<Array<Real>>
    internal var vectorB: Array<Real>
    internal var U: Array<Real>
    internal var I: Array<Real>
    internal var Ignd: Array<Real>

    init {
        m3db = arrayOf<Array<Real>>()
        vectorB = Array(mesh.size()){0.R}
        U = Array(mesh.size()){0.R}
        I = Array(mesh.size()){0.R}
        Ignd = Array(mesh.size()){0.R}
    }
}