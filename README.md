# Drawing Robot

This is my New Providence High School AP Physics C final project.

**Designed to be run in Python 2.7**

### Part 1
The first of this project takes a vector image and returns several sets of coordinates that compose the image. For shapes that are not straight lines, such as cubic and quartic bezier curves, it approximates the curve to within a specified error. In addition, it scales the coordinates to user specified dimensions, calculates boundary coordinates if the paper is too small for the image, and optimizes the ordering off the sets of coordinates to decrease travel distance between disjoint sets of coordinates.

### Part 2
The second part sends instructions to the robot via Bluetooth.

---

**User_Interface.py** - User interface for slicing vector image file into coordinate set

Run with `python User_Interface.py` and follow prompts.

![User Interface](https://github.com/shansteven/drawing-bot/blob/master/Screenshots/example2.png)

<table>
  <tr>
    <td>Original</td>
    <td>Sliced Preview</td>
  </tr>
  <tr>
    <td><img src="https://github.com/shansteven/drawing-bot/blob/master/Screenshots/original.png"></td>
    <td><img src="https://github.com/shansteven/drawing-bot/blob/master/Screenshots/example.png"></td>
  </tr>
</table>

Homer Simpson example from [http://thenewcode.com/assets/images/thumbnails/homer-simpson.svg](http://thenewcode.com/assets/images/thumbnails/homer-simpson.svg)
