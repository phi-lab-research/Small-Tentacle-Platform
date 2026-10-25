# Motion Capture (OptiTrack)

The arm's pose is measured with an OptiTrack system running **Motive** on the Windows
mocap PC. Data is streamed over NatNet to the lab's Ubuntu workstation and published as
ROS 2 topics.

## Main reference

**[tud-phi/ros2-mocap_optitrack](https://github.com/tud-phi/ros2-mocap_optitrack)**:
ROS 2 package for motion capture of soft robots with OptiTrack cameras.

Follow that repository for network setup, build and launch. The rest of this page adds
lab-specific notes on top of it.

## Workflow at a glance

1. Connect the Ubuntu workstation to the mocap network.
2. Open Motive and load the calibration file (or calibrate, the first time).
3. Select the rigid bodies / markers to stream.
4. Launch the ROS 2 package.
5. Read the data from the ROS topics.

---

## 1. Network connection (Ubuntu workstation)

The Motive PC uses the static IP `192.168.4.31` (see the main repo for its setup).

First time only, create a wired connection on Ubuntu:

1. **Settings → Network → Wired → +**
2. **Identity:** choose a name; MTU automatic, leave the rest empty.
3. **IPv4:** method **Manual**
   - Address: `192.168.4.20` (any free last number works)
   - Netmask: `255.255.255.0`
4. **Apply**.

Every time after that, check that this wired connection is the active one.
Test with `ping 192.168.4.31`.

## 2. ROS 2 package

Clone [ros2-mocap_optitrack](https://github.com/tud-phi/ros2-mocap_optitrack) into the
`src/` folder of the ROS 2 workspace (`~/ros2_ws`, see [`code/`](../../code/README.md)) and set the Motive IP in
`mocap_optitrack_client/config/natnetclient.yaml`.

Build and source **every time** you open a terminal:

```bash
cd ~/ros2_ws
colcon build --packages-select mocap_optitrack_client mocap_optitrack_inv_kin mocap_optitrack_w2b
. install/setup.bash
```

Launch from the package folder:

```bash
cd src/ros2-mocap_optitrack
ros2 launch launch/launch_y_up.py
```

> **Note on frames:** in `launch_y_up.py` the `w2b` (world-to-base) node is currently
> commented out. Data is streamed in the Motive world frame, and the transformation to the
> robot base frame is done inside the controller code, so that it stays consistent with the
> soft robot model.

## 3. Motive

### Camera settings

| Setting | Value | Note |
|---|---|---|
| Frame rate | 240 Hz | |
| Exposure | 10 µs | Kept low to reduce reflections |

### Calibration

Usually you do **not** need to recalibrate: open the existing camera calibration file.
Recalibrate only if cameras were moved or tracking is poor.

1. **Mask** existing markers and reflections (*Mask markers*).
2. **Start Wanding** and move the wand through the whole capture volume until every
   camera has at least **3,000 points**.
3. **Calculate**, then **Apply Calibration**.
4. **Clear the masks** (*Clear Mask* button on the right, or per camera: right-click →
   Hardware Mask → Clear). This is needed before the next step, otherwise the ground
   square is not seen.
5. **Ground plane:** place the triangular OptiTrack ground square in the volume. Motive
   recognizes it; select **Auto (passive)** to set the ground plane.
6. **Save** the calibration file.

### Streaming (NatNet)

Follow the main repo, plus:

1. **Edit → Settings → Streaming**, turn streaming **on**.
2. **Local Interface:** `192.168.4.31`; transmission type **Multicast**.
3. Turn on **Labeled markers**, even if you only need unlabeled ones. Without it,
   unlabeled markers are not streamed.

### Assets and rigid bodies

> **Current setup:** only the **world frame** and **rigid objects** are defined as rigid
> bodies. The markers on the soft arm are **not** in a rigid body: they are streamed as
> **free (unlabeled) markers**.
>
> **Check the marker IDs every time you start the controller.** If a marker disappears
> during motion, it comes back with a new, random ID, and the IDs used in the controller
> must be updated.

In **Assets**, select only the assets you want to stream and deselect the others.

**Create a rigid body:**

1. Select the markers → right-click → **Rigid Body → Create From Selected Markers**.
2. Name it in the panel on the right.

> **A marker that belongs to a rigid body no longer appears as a free (unlabeled) marker.**
> Only put markers in a rigid body if you want them tracked as part of that body.

**Rigid-body base frame:** if you use a rigid body as reference frame, give it
**Streaming ID = 1** and make sure the geometric centre of its markers lies where the frame
should be (e.g. on the arm's backbone).

> This applies **only to rigid bodies**. The soft arm's base marker is a **free marker**
> like the other arm markers, so it does not need a rigid body with ID = 1.

**Offsets:** if the body's centre is offset from where it should be, open the
**Builder** pane → select the body → adjust **Orientation** and **Location**.

**Many markers:** if it is hard to place the geometric centre with many markers, create
the rigid body from fewer (e.g. 4) markers first. Then select the body and the remaining
markers → right-click → **Rigid Body → Add Selected Markers to Rigid Body**.

### Recording

Before recording a take, select a working folder with **New → Open Folder**.

---

## Things to know

- **Rigid bodies need rigid markers:** if one marker in a rigid body moves relative to the
  others (e.g. on the soft arm), the body is lost.
- **Body frames** are shown only for rigid bodies.
- **Marker-only tracking in Motive:** if the body's markers are activated in Motive, the
  body is shown as still even when the markers move.
- **Manual labels and marker sets** can only be assigned to recordings, not live.
- **Unlabeled marker IDs are not stable:** if a marker is lost from view its ID changes,
  and the order changes with every recording.

## Markers on the arm

Marker positions along the arm are in
[`data/soft_arm_parameters.yaml`](../../data/soft_arm_parameters.yaml) under `mocap_markers`.
