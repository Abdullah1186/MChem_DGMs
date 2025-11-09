from flask import Flask, render_template

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('index.html')  # looks in templates/index.html


@app.route('/analysis')
def analysis():
    return render_template('analysis.html') # going here so u can run analysis e.g return structural_analysis()


if __name__ == "__main__":
    app.run(debug=True)