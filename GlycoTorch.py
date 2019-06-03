from flask import Flask
from flask import render_template
from flask import flash, request, redirect
from werkzeug.utils import secure_filename
from flask import send_from_directory
from flask_cors import CORS, cross_origin
from Carbohydrate import *
from Docking_Analysis import *

app = Flask(__name__)


UPLOAD_FOLDER = '/tmp'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# DATA_FOLDER = '/data/ligands/pdb'
# app.config['DATA_FOLDER'] = DATA_FOLDER

# GLYCOSIDIC_FOLDER = '/notebooks/glycosidic'
# app.config['GLYCOSIDIC_FOLDER'] = GLYCOSIDIC_FOLDER

app.config['DEBUG'] = True

CORS(app)

ALLOWED_EXTENSIONS = {'pdb', 'pdbqt'}


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/', methods=['GET', 'POST'])
def upload_file():

    if request.method == 'POST':
        print(request.form.__dict__)
        if 'analyse_ligand' in request.form:

            # check if the post request has the file part
            if 'ligand_to_analyse' not in request.files:
                flash('No file part')
                return redirect(request.url)

            file = request.files['ligand_to_analyse']
            # if user does not select file, browser also
            # submit an empty part without filename
            if file.filename == '':
                flash('No selected file')
                return redirect(request.url)
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
                return redirect('/LigandAnalysis/'+filename)

        if 'docking_analysis' in request.form:
            # check if the post request has the file part
            if 'vina_output_file' not in request.files:
                flash('No vina_output_file found')
                return redirect(request.url)
            if 'vina_protein_file' not in request.files:
                flash('No vina_protein_file found')
                return redirect(request.url)

            vina_output_file = request.files['vina_output_file']
            vina_protein_file = request.files['vina_protein_file']
            vina_crystal_pose = request.files['vina_crystal_pose']

            # if user does not select file, browser also
            # submit an empty part without filename
            if file.filename == '':
                flash('No selected file')
                return redirect(request.url)

            if vina_output_file \
                    and allowed_file(vina_output_file.filename) and vina_protein_file and \
                    allowed_file(vina_protein_file.filename) and vina_crystal_pose and \
                    allowed_file(vina_crystal_pose.filename):


                vina_output_filename = secure_filename(vina_output_file.filename)
                vina_protein_filename = secure_filename(vina_protein_file.filename)
                vina_crystal_pose_filename = secure_filename(vina_crystal_pose.filename)


                file.save(os.path.join(app.config['UPLOAD_FOLDER'], vina_output_filename))
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], vina_protein_filename))
                file.save(os.path.join(app.config['UPLOAD_FOLDER'], vina_crystal_pose_filename))

                tmp = "/home/EricBoittier/mysite/"

                docking_analysis(tmp+"/"+vina_protein_filename, tmp+"/"+vina_output_filename, tmp+"/"+vina_crystal_pose_filename, location=tmp)


                return render_template('Inbound_Docking_Analysis.html', )

    return render_template('index.html')

@app.route('/tmp/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'],
                               filename)

@app.route("/summary_of_results/")
def summary_of_results():
    l = os.listdir("/home/EricBoittier/mysite/data/results/glycosidic/0_0/small")
    return render_template("results_index.html", list=l)


@app.route("/vinacarbresult/<path1>/<path2>/<path3>")
def vinacarb(path1, path2, path3):
    path = path1 + "/" + path2 + "/" + path3

    cluster_text_files = os.listdir("/home/EricBoittier/mysite/data/results/glycosidic/"+path)
    cluster_text_files = [c for c in cluster_text_files if c.__contains__("cluster_")]

    crystal = path3 + ".pdb.mol2.pdbqt.pdb"

    protein = os.listdir("/home/EricBoittier/mysite/data/results/glycosidic/"+path)
    protein = [c for c in protein if c.__contains__("WithoutLigands.pdb.pdbqt.pdb")]

    clusters = {}
    for text_file in cluster_text_files:
        temp = []
        f = open("/home/EricBoittier/mysite/data/results/glycosidic/"+path+"/"+text_file, "r")
        for l in f.readlines():
            temp.append(l)
        clusters[text_file] = temp

    barchart = [b for b in os.listdir("/home/EricBoittier/mysite/data/results/glycosidic/"+path) if b.__contains__("_barchart_")]


    return render_template("Local_Docking_Analysis.html", clusters=clusters,
                           barchart=barchart[0],
                           crystal_pose=crystal, protein=protein[0], path=path)

@app.route("/results/<type>")
def get_results(type):
    pdfs = os.listdir( "/home/EricBoittier/mysite/notebooks/" + type )
    return render_template("results.html", type=type, pdfs=pdfs)

@app.route('/data/ligands/pdb/<filename>')
def get_ligand_pdb(filename):
    return send_from_directory(app.config['DATA_FOLDER'],
                               filename)

@app.route('/notebooks/glycosidic/<filename>')
def get_glycosidic(filename):
    return send_from_directory(app.config['GLYCOSIDIC_FOLDER'],
                               filename)

@app.route('/uploads/uploads.html')
def get_uploads():
    return render_template('/uploads/uploads.html')


@app.route('/LigandAnalysis/<name>')
def ligandAnalysis(name):
    filename = os.path.join(app.config['UPLOAD_FOLDER'], name)
    l = Carbohydrate(filename)

    if request.method == 'POST':
        if 'download' in request.form::
            # hello
            flash('Working on it')
    elif request.method == 'GET':
        return render_template('LigandAnalysis.html', name=name, ligand=l)


